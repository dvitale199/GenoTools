# From: https://github.com/GenoML/genoml2
#
# Copyright 2020 The GenoML Authors. All Rights Reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# ==============================================================================

import io
import logging
import os
import pathlib
import platform
import requests
import stat
import subprocess
import zipfile
import tarfile
import sys

# can decide where to throw the executable files later
## for now /$HOME/.genotools/misc/executables
def __get_executable_folder():
    key = "GENOTOOLS_DEP_DIR"
    if key in os.environ:
        return os.path.abspath(os.environ.get(key))
    else:
        return os.path.join(str(pathlib.Path.home()), ".genotools", "misc",
                            "executables")


__executable_folder = __get_executable_folder()

def shell_do(command, log=False, return_log=False):
    print(f'Executing: {(" ").join(command.split())}', file=sys.stderr)

    res=subprocess.run(command.split(), stdout=subprocess.PIPE)

    if log:
        print(res.stdout.decode('utf-8'))
    if return_log:
        return(res.stdout.decode('utf-8'))

def __check_exec(exec_path, *args, absolute_path=False):
    if not absolute_path:
        binary_path = os.path.join(__executable_folder, exec_path)
    else:
        binary_path = exec_path
    if not os.path.exists(binary_path):
        return False

    _ = subprocess.run([binary_path, *args], stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL)
    return True


def __install_exec(url, exec_path):
    r = requests.get(url, verify=False, stream=True)

    if '.zip' in url:
        r.raw.decode_content = True
        buffer = io.BytesIO()
        buffer.write(r.content)
        with zipfile.ZipFile(buffer, "r") as fp:
            fp.extractall(__executable_folder)
        
    elif '.tar.gz' in url:
        file = tarfile.open(fileobj=r.raw, mode="r|gz")
        file.extractall(__executable_folder)

    binary_path = os.path.join(__executable_folder, exec_path)

    os.chmod(binary_path, stat.S_IEXEC)


def __check_package(name):
    platform_system = platform.system()

    if name not in __DEPENDENCIES:
        raise EnvironmentError("Unknown package: {}".format(name))

    if platform_system not in __DEPENDENCIES[name]:
        raise EnvironmentError(
            "Unknown supported OK: {}".format(platform_system))

    entry = __DEPENDENCIES[name][platform_system]

    binary_name = entry["binary"]
    args = entry["version_args"]
    url = entry["url"]

    if __check_exec(binary_name, *args):
        print(f'{name} is found')
        logging.debug("{} is found".format(name))
        return os.path.join(__executable_folder, binary_name)

    logging.warning("Installing {}".format(name))
    __install_exec(url, binary_name)
    if not __check_exec(binary_name, *args):
        logging.warning("Failed to run {} after installation".format(name))
        raise EnvironmentError("Can not install {}".format(name))
    else:
        return os.path.join(__executable_folder, binary_name)


def check_dependencies():
    global __DEPENDENCIES
    ret = {}
    for package, data in __DEPENDENCIES.items():
        if "checker" in data:
            # instead of GenoMLs description loader
            print('Loaded')

    return ret

def check_plink():
    return __check_package('Plink')

def check_plink2():
    return __check_package('Plink2')

__DEPENDENCIES = {
    'Plink': {
        'checker': check_plink,
        'Darwin': {
            'binary': 'plink',
            'version_args': ['--version'],
            'url': 'https://s3.amazonaws.com/plink1-assets/dev/plink_mac.zip'
        },
        'Linux': {
            'binary': 'plink',
            'version_args': ['--version'],
            'url': 'https://s3.amazonaws.com/plink1-assets/dev/plink_linux_x86_64.zip'
        },
        'Windows': {
            'binary': 'plink',
            'version_args': ['--version'],
            'url': 'https://s3.amazonaws.com/plink1-assets/plink_win64_20220402.zip'
        }
    },

    'Plink2': {
        'checker': check_plink2,
        'Darwin': {
            'binary': 'plink2',
            'version_args': ['--version'],
            'url': 'https://s3.amazonaws.com/plink2-assets/alpha3/plink2_mac_avx2_20220603.zip'
        },
        'Linux': {
            'binary': 'plink2',
            'version_args': ['--version'],
            'url': 'https://s3.amazonaws.com/plink2-assets/alpha3/plink2_linux_x86_64_20220603.zip'
        }
    }
}
