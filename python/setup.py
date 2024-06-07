# Md Fahad Hasan, 2023, 2024
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; version 3, or (at your option) any
# later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
#


from setuptools import setup, Extension, find_packages, Distribution
import os
import platform
import importlib.metadata
import sys
import subprocess
import filecmp
from distutils.cmd import Command

# Boost library version based on Python version
python_major_version = sys.version_info.major
python_minor_version = sys.version_info.minor
boost_python_lib = f"boost_python{python_major_version}{python_minor_version}"
boost_numpy_lib = f"boost_numpy{python_major_version}{python_minor_version}"

class CreateSymlinks(Command):
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        target_dirs = {
            'src': 'src',
            'treedec': 'treedec',
            'tests': 'tests',
            'tdlib': 'tdlib',
            'examples': 'examples',
        }
        exclude_file = '__init__.py'
        header_files = [
            'container.hpp',
            'graph_impl.hpp',
            'graph_traits.hpp',
            'treedec_traits.hpp',
        ]

        for target, link_name in target_dirs.items():
            target_dir = os.path.join(root_dir, target)
            link_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), link_name)
            if not os.path.exists(link_dir):
                if target == 'tdlib':
                    os.makedirs(link_dir, exist_ok=True)
                    for item in os.listdir(target_dir):
                        if item != exclude_file:
                            source_path = os.path.join(target_dir, item)
                            link_path = os.path.join(link_dir, item)
                            if not os.path.exists(link_path):
                                os.symlink(source_path, link_path)
                else:
                    os.symlink(target_dir, link_dir)

            if target == 'treedec':
                treedec_subfolder = os.path.join(link_dir, 'treedec')
                if not os.path.exists(treedec_subfolder):
                    os.makedirs(treedec_subfolder, exist_ok=True)

                src_dir = os.path.join(root_dir, 'src')
                for header_file in header_files:
                    src_file_path = os.path.join(src_dir, header_file)
                    symlink_path = os.path.join(treedec_subfolder, header_file)
                    if not os.path.exists(symlink_path):
                        os.symlink(src_file_path, symlink_path)


def find_include_path_in_package(package_name):
    try:
        distribution_files = importlib.metadata.files(package_name)
        record_file = next(f for f in distribution_files if f.name == 'RECORD')
        record_file_path = record_file.locate()

        with open(record_file_path, 'r') as file:
            record_contents = file.read()

        for line in record_contents.splitlines():
            path, _, _ = line.partition(',')
            if 'include' in path and path.endswith('boost.h'):
                base_path = os.path.dirname(record_file_path)
                full_path = os.path.join(base_path, '..', path)
                normalized_path = os.path.normpath(full_path)
                absolute_path = os.path.abspath(normalized_path)
                if os.path.exists(absolute_path):
                    return absolute_path.replace('/gala/boost.h', '')

    except StopIteration:
        pass

    return None

include_path = find_include_path_in_package('gala')
include_dirs = ['treedec', "src"]
if include_path:
    print(f"Gala path found at: {include_path}")
    include_dirs.append(include_path)
else:
    print("Include path not found for package 'gala'.")


def get_boost_python_dirs_for_macos():
    try:
        brew_prefix = subprocess.check_output(['brew', '--prefix']).strip().decode()

        boost_python_cellar_path = os.path.join(brew_prefix, 'Cellar', 'boost-python3')
        boost_python_versions = [d for d in os.listdir(boost_python_cellar_path) if os.path.isdir(os.path.join(boost_python_cellar_path, d))]
        if not boost_python_versions:
            return None, None

        latest_version = max(boost_python_versions)
        lib_dir = os.path.join(boost_python_cellar_path, latest_version, 'lib')

        boost_include_path = os.path.join(brew_prefix, 'include')
        
        return lib_dir, boost_include_path
    except Exception as e:
        print(f"An error occurred: {e}")
        return None, None


boost_python_lib_dir_for_macos, boost_python_include_dir_for_macos = None, None
if platform.system() == 'Darwin':
    boost_python_lib_dir_for_macos, boost_python_include_dir_for_macos = get_boost_python_dirs_for_macos()
    if boost_python_lib_dir_for_macos and boost_python_include_dir_for_macos:
        include_dirs.append(boost_python_include_dir_for_macos)
        print(f"Boost-Python library directory: {boost_python_lib_dir_for_macos}")
        print(f"Boost include directory: {boost_python_include_dir_for_macos}")
    else:
        print("Boost-Python not found via Homebrew")

def fetch_from_makelist(file_path):
    with open(file_path, 'r') as file:
        # BUG: read line by line.
        content = file.read()

    libraries = []
    modules = []

    # BUG: sensitive to whitespace etc.
    for line in content.splitlines():
        if line.endswith('.la \\'):
            libraries.append("treedec/" + line[:-5].strip())
            modules.append(line[:-5].strip())
        elif line.endswith('.la'):
            libraries.append("treedec/" + line[:-3].strip())
            modules.append(line[:-3].strip())
        
    return libraries, modules

makelist_path = '../treedec/MakeList'
cpp_la, modules = fetch_from_makelist(makelist_path)
cpp_src = [i + ".cpp" for i in cpp_la]

    
def check_lib_folder(build_dir):
    if not os.path.exists(build_dir):
        print(f"Build directory '{build_dir}' does not exist.")
        return
    files = os.listdir(build_dir)
    lib_folder = None
    for file in files:
        if file.startswith('lib.'):
            lib_folder = file
            break
    if lib_folder:
        print(f"Found lib folder '{lib_folder}' in the build directory.")
        return lib_folder
    else:
        print("No lib folder found in the build directory.")
    return None

def create_test_command_class(test_directories, reference_dir):
    class TestCommand(Command):
        user_options = []

        def initialize_options(self):
            pass

        def finalize_options(self):
            pass

        def run(self):
            lib_folder = check_lib_folder('build')
            os.environ['PYTHONPATH'] = '.' + 'treedec/tests:treedec/tests/tdlib:treedec/tests:treedec/tests/tdlib:tdlib:tdlib/tests:build/'+str(lib_folder)
            output_dir = "test_outputs"
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
            
            for directory_prefix, test_files in test_directories.items():
                for test_file in test_files:
                    full_test_file_path = os.path.join(directory_prefix, test_file)
                    test_name = os.path.basename(test_file)
                    output_file_path = os.path.join(output_dir, f"{test_name}.out")
                    reference_file_name = f"{test_name}.out"
                    reference_file_path = os.path.join(reference_dir, reference_file_name)

                    if test_file == "test_long_approx_CFGs.py":
                        command = f"python{python_major_version}.{python_minor_version} {full_test_file_path} long"
                    else:
                        command = f"python{python_major_version}.{python_minor_version} {full_test_file_path}"

                    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=sys.stderr)
                    stdout, stderr = process.communicate()

                    with open(output_file_path, 'w') as f:
                        f.write(stdout.decode())
            
                    try:
                        if (process.returncode != 0 and process.returncode != 1) or not os.path.exists(reference_file_path) or not filecmp.cmp(output_file_path, reference_file_path):
                            print(f"Test {full_test_file_path} failed")
                            if process.returncode == 77:
                                print("Warning: Need to fix")
                            if process.returncode != 0:
                                print("Return code:", process.returncode)
                            if not filecmp.cmp(output_file_path, reference_file_path, shallow=False):
                                print("Output does not match reference")
                        else:
                            print(f"Test {full_test_file_path} passed")
                    except FileNotFoundError:
                        print(f"Reference file not found: {reference_file_path}, test {full_test_file_path} cannot be verified")
    return TestCommand

test_directories = {
    "../tdlib/tests": [
        "test_pp_fi.py", "exact_ta.py", "dclc_g5.py", "test_long_fi_Dimacs.py", "test_long_approx_CFGs.py"
    ],
    "../treedec/tests": [
        "graph_bal0.py", "graph_gala0.py", "excut_0.py", "exta_0.py", "fi_0.py", "fi_1.py", "fi_dimacs.py", "order0.py", "pp_dimacs.py", "pp_cfg.py", "misc_0.py"
    ]
}

TestCommand = create_test_command_class(test_directories, "treedec/tests/ref")

distribution = Distribution()
create_symlinks_cmd = CreateSymlinks(distribution)
create_symlinks_cmd.initialize_options()
create_symlinks_cmd.finalize_options()
create_symlinks_cmd.run()

extensions = [
    Extension(
        'treedec.{}'.format(modules[i]), 
        [cpp_src[i]],
        include_dirs=include_dirs,
        language='c++',
        libraries=[boost_python_lib, boost_numpy_lib],
        extra_compile_args=['-std=c++20', '-fmax-errors=1', '-DHAVE_GALA_GRAPH_H'],
        **({'library_dirs': [boost_python_lib_dir_for_macos]} if platform.system() == 'Darwin' and boost_python_lib_dir_for_macos != None else {}) 
    )
    for i, source in enumerate(cpp_src)
]

setup(
    name='treedec',
    version='0.9.2+dev',
    packages=find_packages(),
    ext_modules=extensions,
    install_requires=['gala'],
    cmdclass={
        'create_symlinks': CreateSymlinks,
        'test': TestCommand,
    }
)
