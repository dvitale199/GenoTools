"""Immutable genotype data representation.

This module provides the GenotypeData class, an immutable reference to
genotype files that handles format detection, conversion, and counting.

Key features:
- Immutable (frozen dataclass) - transformations return new instances
- Format-aware - detects pfile vs bfile automatically
- Lazy conversion - only converts format when necessary
- Consistent path handling - paths never include extensions

Usage:
    from genotools.core.genotypes import GenotypeData

    # Load from existing files (auto-detects format)
    data = GenotypeData.from_path(Path("/path/to/genotypes"))

    # Convert to different format (returns NEW instance)
    pfile_data = data.to_pfile(Path("/path/to/output"))

    # Ensure format for operations that require it
    bfile_data = data.ensure_bfile(work_dir)
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Literal, Optional, Tuple
import subprocess

from .exceptions import FileFormatError, ValidationError, ExternalToolError
from .logging import get_logger

logger = get_logger(__name__)

# File extensions for each format
PFILE_EXTENSIONS = (".pgen", ".pvar", ".psam")
BFILE_EXTENSIONS = (".bed", ".bim", ".fam")

# CRITICAL: Always use this flag to preserve sample metadata in psam files
PLINK2_PSAM_COLS = "psam-cols=fid,parents,sex,pheno1,phenos"


@dataclass(frozen=True)
class GenotypeData:
    """Immutable reference to genotype files.

    This class represents a set of genotype files (either pfile or bfile format)
    and provides methods for format detection, conversion, and basic statistics.

    Attributes:
        path: Path to the genotype files (without extension)
        format: Either "pfile" (PLINK2) or "bfile" (PLINK1.9)
        sample_count: Number of samples in the dataset
        variant_count: Number of variants in the dataset

    Note:
        This is a frozen dataclass - all attributes are immutable.
        Transformation methods (to_pfile, to_bfile) return NEW instances.
    """

    path: Path
    format: Literal["pfile", "bfile"]
    sample_count: int
    variant_count: int

    @classmethod
    def from_path(cls, path: Path | str) -> "GenotypeData":
        """Create GenotypeData from existing files.

        Automatically detects the format (pfile or bfile) and counts
        samples and variants.

        Args:
            path: Path to genotype files (without extension).
                  Can be string or Path object.

        Returns:
            GenotypeData instance with detected format and counts.

        Raises:
            FileFormatError: If neither pfile nor bfile format is found.
            FileNotFoundError: If required files are missing.

        Example:
            data = GenotypeData.from_path("/data/my_genotypes")
            # Looks for my_genotypes.pgen/.pvar/.psam or .bed/.bim/.fam
        """
        path = Path(path)

        # Detect format
        format_type = cls._detect_format(path)

        # Count samples and variants
        sample_count, variant_count = cls._count_samples_variants(path, format_type)

        logger.debug(
            f"Loaded GenotypeData: {path}, format={format_type}, "
            f"samples={sample_count}, variants={variant_count}"
        )

        return cls(
            path=path,
            format=format_type,
            sample_count=sample_count,
            variant_count=variant_count
        )

    @staticmethod
    def _detect_format(path: Path) -> Literal["pfile", "bfile"]:
        """Detect genotype file format from existing files.

        Args:
            path: Base path (without extension)

        Returns:
            "pfile" or "bfile"

        Raises:
            FileFormatError: If neither format is found
        """
        # Check for pfile format first (preferred)
        pfile_exists = all(
            path.with_suffix(ext).exists() for ext in PFILE_EXTENSIONS
        )
        if pfile_exists:
            return "pfile"

        # Check for bfile format
        bfile_exists = all(
            path.with_suffix(ext).exists() for ext in BFILE_EXTENSIONS
        )
        if bfile_exists:
            return "bfile"

        # Neither format found
        raise FileFormatError(
            f"No valid genotype files found at {path}. "
            f"Expected either pfile (.pgen/.pvar/.psam) or bfile (.bed/.bim/.fam) format."
        )

    @staticmethod
    def _count_samples_variants(
        path: Path,
        format_type: Literal["pfile", "bfile"]
    ) -> Tuple[int, int]:
        """Count samples and variants in genotype files.

        Args:
            path: Base path (without extension)
            format_type: Either "pfile" or "bfile"

        Returns:
            Tuple of (sample_count, variant_count)
        """
        if format_type == "pfile":
            sample_file = path.with_suffix(".psam")
            variant_file = path.with_suffix(".pvar")

            # Count samples (subtract 1 for header)
            sample_count = sum(1 for _ in open(sample_file)) - 1

            # Count variants (handle comment lines in pvar)
            variant_count = 0
            with open(variant_file) as f:
                for line in f:
                    if not line.startswith("#"):
                        variant_count += 1

        else:  # bfile
            sample_file = path.with_suffix(".fam")
            variant_file = path.with_suffix(".bim")

            # fam has no header
            sample_count = sum(1 for _ in open(sample_file))
            # bim has no header
            variant_count = sum(1 for _ in open(variant_file))

        return sample_count, variant_count

    def to_pfile(
        self,
        out_path: Path | str,
        plink2_path: Optional[Path | str] = None
    ) -> "GenotypeData":
        """Convert to pfile format.

        If already in pfile format and out_path is different from current path,
        creates a copy. If already pfile and same path, returns self.

        Args:
            out_path: Output path (without extension)
            plink2_path: Path to plink2 executable. If None, will be found
                        automatically via dependencies module.

        Returns:
            NEW GenotypeData instance pointing to the pfile output.

        Raises:
            ExternalToolError: If PLINK2 conversion fails.
        """
        out_path = Path(out_path)

        # If already pfile and same path, return self (immutable)
        if self.format == "pfile" and out_path == self.path:
            return self

        # Get plink2 executable
        if plink2_path is None:
            from genotools.dependencies import check_plink2
            plink2_path = check_plink2()  # type: ignore[no-untyped-call]

        # Build conversion command
        if self.format == "pfile":
            input_flag = "--pfile"
        else:
            input_flag = "--bfile"

        # Ensure output directory exists
        out_path.parent.mkdir(parents=True, exist_ok=True)

        cmd = [
            str(plink2_path),
            input_flag, str(self.path),
            "--make-pgen", PLINK2_PSAM_COLS,
            "--out", str(out_path)
        ]

        logger.debug(f"Converting to pfile: {' '.join(cmd)}")

        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            raise ExternalToolError(
                tool="plink2",
                command=" ".join(cmd),
                returncode=result.returncode,
                stderr=result.stderr,
                stdout=result.stdout
            )

        # Clean up log file (following existing pattern)
        log_file = out_path.with_suffix(".log")
        if log_file.exists():
            log_file.unlink()

        return GenotypeData.from_path(out_path)

    def to_bfile(
        self,
        out_path: Path | str,
        plink2_path: Optional[Path | str] = None
    ) -> "GenotypeData":
        """Convert to bfile format.

        If already in bfile format and out_path is different from current path,
        creates a copy. If already bfile and same path, returns self.

        Args:
            out_path: Output path (without extension)
            plink2_path: Path to plink2 executable. If None, will be found
                        automatically via dependencies module.

        Returns:
            NEW GenotypeData instance pointing to the bfile output.

        Raises:
            ExternalToolError: If PLINK2 conversion fails.
        """
        out_path = Path(out_path)

        # If already bfile and same path, return self (immutable)
        if self.format == "bfile" and out_path == self.path:
            return self

        # Get plink2 executable
        if plink2_path is None:
            from genotools.dependencies import check_plink2
            plink2_path = check_plink2()  # type: ignore[no-untyped-call]

        # Build conversion command
        if self.format == "pfile":
            input_flag = "--pfile"
        else:
            input_flag = "--bfile"

        # Ensure output directory exists
        out_path.parent.mkdir(parents=True, exist_ok=True)

        cmd = [
            str(plink2_path),
            input_flag, str(self.path),
            "--make-bed",
            "--out", str(out_path)
        ]

        logger.debug(f"Converting to bfile: {' '.join(cmd)}")

        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            raise ExternalToolError(
                tool="plink2",
                command=" ".join(cmd),
                returncode=result.returncode,
                stderr=result.stderr,
                stdout=result.stdout
            )

        # Clean up log file (following existing pattern)
        log_file = out_path.with_suffix(".log")
        if log_file.exists():
            log_file.unlink()

        return GenotypeData.from_path(out_path)

    def ensure_pfile(self, work_dir: Path | str) -> "GenotypeData":
        """Ensure data is in pfile format, converting if necessary.

        If already pfile, returns self. Otherwise converts to pfile
        in the specified work directory.

        Args:
            work_dir: Directory for conversion output if needed.

        Returns:
            GenotypeData in pfile format (self if already pfile).
        """
        if self.format == "pfile":
            return self

        work_dir = Path(work_dir)
        out_path = work_dir / f"{self.path.stem}_pfile"
        return self.to_pfile(out_path)

    def ensure_bfile(self, work_dir: Path | str) -> "GenotypeData":
        """Ensure data is in bfile format, converting if necessary.

        If already bfile, returns self. Otherwise converts to bfile
        in the specified work directory.

        Args:
            work_dir: Directory for conversion output if needed.

        Returns:
            GenotypeData in bfile format (self if already bfile).
        """
        if self.format == "bfile":
            return self

        work_dir = Path(work_dir)
        out_path = work_dir / f"{self.path.stem}_bfile"
        return self.to_bfile(out_path)

    def validate_required_columns(self) -> None:
        """Validate that required columns exist in sample file.

        Checks for SEX and PHENO1 columns which are required for
        many QC operations.

        Raises:
            ValidationError: If required columns are missing.
        """
        if self.format == "pfile":
            sample_file = self.path.with_suffix(".psam")
            with open(sample_file) as f:
                header = f.readline().strip()

            # psam header starts with #FID or #IID
            columns = header.replace("#", "").split()

            if "SEX" not in columns:
                raise ValidationError(
                    f"{sample_file} is missing SEX column. "
                    "GenoTools requires a SEX column even if no sex data is available."
                )

            if "PHENO1" not in columns:
                raise ValidationError(
                    f"{sample_file} is missing PHENO1 column. "
                    "GenoTools requires a PHENO1 column even if no phenotype data is available."
                )
        else:
            # bfile format - fam file has fixed columns, no header
            # Columns: FID IID PAT MAT SEX PHENO
            # These are always present by format definition
            pass

    @property
    def pgen(self) -> Optional[Path]:
        """Path to .pgen file if pfile format, else None."""
        if self.format == "pfile":
            return self.path.with_suffix(".pgen")
        return None

    @property
    def pvar(self) -> Optional[Path]:
        """Path to .pvar file if pfile format, else None."""
        if self.format == "pfile":
            return self.path.with_suffix(".pvar")
        return None

    @property
    def psam(self) -> Optional[Path]:
        """Path to .psam file if pfile format, else None."""
        if self.format == "pfile":
            return self.path.with_suffix(".psam")
        return None

    @property
    def bed(self) -> Optional[Path]:
        """Path to .bed file if bfile format, else None."""
        if self.format == "bfile":
            return self.path.with_suffix(".bed")
        return None

    @property
    def bim(self) -> Optional[Path]:
        """Path to .bim file if bfile format, else None."""
        if self.format == "bfile":
            return self.path.with_suffix(".bim")
        return None

    @property
    def fam(self) -> Optional[Path]:
        """Path to .fam file if bfile format, else None."""
        if self.format == "bfile":
            return self.path.with_suffix(".fam")
        return None

    def __str__(self) -> str:
        return (
            f"GenotypeData({self.path}, format={self.format}, "
            f"samples={self.sample_count}, variants={self.variant_count})"
        )
