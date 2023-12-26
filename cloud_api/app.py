import streamlit as st
import subprocess

# Header
st.title("Genotools QC and Ancestry Analysis")

# Dropdown for model selection
selected_model = st.selectbox("Select model", ["model1", "model2", "model3"])

# Input paths
pfile_path = st.text_input("Path to genotypes for QC")
out_path = st.text_input("Path to QC output")
ref_panel_path = st.text_input("Path to reference panel")
ref_labels_path = st.text_input("Path to reference ancestry labels")

# Run button
if st.button("Run Genotools"):
    # Construct the command with model option
    command = f"genotools --pfile {pfile_path} --out {out_path} --ancestry --ref_panel {ref_panel_path} --ref_labels {ref_labels_path} --all_sample --all_variant --model /path/to/{selected_model}"

    # Run the command and display output
    try:
        result = subprocess.run(command, shell=True, check=True, capture_output=True, text=True)
        st.write("Genotools command executed successfully!")
        st.write(result.stdout)  # Display command output
    except subprocess.CalledProcessError as e:
        st.error("Genotools command failed with error:")
        st.error(e.output)  # Display error output