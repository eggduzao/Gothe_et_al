import os
import sys
import io
import re
from shutil import copy
from pwd import getpwnam
from sys import platform, exit
from distutils import dir_util
from setuptools import setup, find_packages
from os import walk, chown, chmod, path, getenv, makedirs, remove
from optparse import OptionParser, BadOptionError, AmbiguousOptionError

# Python 3 support
#if not sys.version_info[0] == 2:
#    sys.exit("Sorry, Python 3 is not supported (yet)")

"""
Installs the Gothe et al. tools with standard setuptools options.

Authors: Eduardo G. Gusmao, Giuseppe Petrosino, Argyris Papantonis, Vassilis Roukos.

Installs the Gothe et al. tools with standard setuptools options.
"""

def read(*names, **kwargs):
    with io.open(
        os.path.join(os.path.dirname(__file__), *names),
        encoding=kwargs.get("encoding", "utf8")
    ) as fp:
        return fp.read()


def find_version(*file_paths):
    version_file = read(*file_paths)
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]",
                              version_file, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")

current_version = find_version("src", "__version__.py")

###################################################################################################
# Unsupported Platforms
###################################################################################################

# Windows is not supported
supported_platforms = ["linux", "linux2", "darwin"]
if platform not in supported_platforms:
    print("ERROR: This package currently supports only unix-based systems (Linux and MAC OS X).")
    exit(0)

###################################################################################################
# Parameters
###################################################################################################

"""
Tools Dictionary:
  * Insert in the dictionary bellow a key+tuple X: (Y,Z,W,K) representing:
    - X: A string representing the name the user should provide for this script in order to install the tool.
    - Y: A string representing the name of the program which the user should type in the terminal in order to execute the tool.
    - Z: A string representing the path to the main function that executes the program.
    - W: A list with the package requirements for that program.
    - K: A list with external binary files that will be copied by the setup installation function to the user's bin folder.
"""

common_deps = ["numpy>=1.4.0",
               "scipy>=0.7.0",
               "pysam>=0.8.2",
               "pyBigWig"]

if platform.startswith("darwin"):
    bin_dir = "mac"
    libG = "lib_mac.so"
else:
    bin_dir = "linux"
    libG = "lib_linux.so"

tools_dictionary = {
"correlation_bliss_features": (
    "corr-dsb-feat",
    "src.correlation_bliss_features.Main:main",
    [],
    []
),
"correlation_dsb_distance_expression": (
    "corr-dsb-dist-exp",
    "src.correlation_dsb_distance_expression.Main:main",
    [],
    []
),
"ctcf_loop_analysis": (
    "ctcf-loop",
    "src.ctcf_loop_analysis.Main:main",
    [],
    []
),
"gene_metaplots": (
    "gene-meta",
    "src.gene_metaplots.Main:main",
    [],
    []
),
"heatmaps_bliss_features": (
    "heatmap-dsb",
    "src.heatmaps_bliss_features.Main:main",
    [],
    []
)
}

###################################################################################################
# Auxiliary Functions/Classes
###################################################################################################

# PassThroughOptionParser Class
class PassThroughOptionParser(OptionParser):
    """
    An unknown option pass-through implementation of OptionParser.
    When unknown arguments are encountered, bundle with largs and try again,
    until rargs is depleted.
    sys.exit(status) will still be called if a known argument is passed
    incorrectly (e.g. missing arguments or bad argument types, etc.)
    """
    def _process_args(self, largs, rargs, values):
        while rargs:
            try:
                OptionParser._process_args(self, largs, rargs, values)
            except (BadOptionError,AmbiguousOptionError), e:
                largs.append(e.opt_str)


# recursive_chown_chmod Function
def recursive_chown_chmod(path_to_walk, uid, gid, file_permission, path_permission):
    """
    Recursively applies chown from path.
    """
    for root_dir, directory_list, file_list in walk(path_to_walk):
        chown(root_dir, uid, gid)
        chmod(root_dir, path_permission)
        for f in file_list:
            current_complete_file = path.join(root_dir,f)
            chown(current_complete_file, uid, gid)
            chmod(current_complete_file, file_permission)

###################################################################################################
# Processing Input Arguments
###################################################################################################

# Parameters
data_base_name = "gothe_et_al"
usage_message = "python setup.py install [python options] [gothe_et_al options]"
version_message = "Spatial chromosome folding and active transcription drive DNA fragility and formation of oncogenic MLL translocations. Version: "+str(current_version)

# Initializing Option Parser
parser = PassThroughOptionParser(usage=usage_message, version=version_message)

# Parameter: Data Location
param_data_location_name = "--local-data-path"
parser.add_option(param_data_location_name, type="string", metavar="STRING",
                  help="Path containing data used by Gothe et al. tool.",
                  dest="param_data_location", default=path.join(getenv('HOME'),data_base_name))

# Parameter: Tool
param_tools_name = "--instal-tools"
parser.add_option(param_tools_name, type="string", metavar="STRING",
                  help=("The tools which will be installed. If this argument is not used, "
                        "then all command-line tools are installed. The current available options "
                        "are: "+", ".join(tools_dictionary.keys())+"; You can also provide " 
                        "multiple tools in a list separated by comma."),
                  dest="param_tools", default=",".join(tools_dictionary.keys()))

# Processing Options
options, arguments = parser.parse_args()
if path.basename(options.param_data_location) != data_base_name:
    options.param_data_location = path.join(options.param_data_location,data_base_name)
if options.param_data_location[0] == "~":
    options.param_data_location = path.join(getenv('HOME'),options.param_data_location[2:])
options.param_tools = options.param_tools.split(",")

# Manually Removing Additional Options from sys.argv
new_sys_argv = []
for e in sys.argv:
    if param_data_location_name == e[:len(param_data_location_name)]:
        continue
    elif param_tools_name == e[:len(param_tools_name)]:
        continue
    new_sys_argv.append(e)

sys.argv = new_sys_argv

# Defining entry points
current_entry_points = {"console_scripts" : []}
for tool_option in options.param_tools:
    if tool_option != "core":
        current_entry_points["console_scripts"].append(" = ".join(tools_dictionary[tool_option][:2]))

# Defining install requirements
current_install_requires = common_deps
for tool_option in options.param_tools:
    current_install_requires += tools_dictionary[tool_option][2]

###################################################################################################
# Setup Function
###################################################################################################

# Parameters
short_description = "Toolkit to perform Gothe et al paper's analyses"
classifiers_list = ["Topic :: Scientific/Engineering :: Bio-Informatics",
                    "Topic :: Scientific/Engineering :: Artificial Intelligence"]
keywords_list = ["ChIP-seq", "DNase-seq", "DNA Damage", "DSBs", "CTCF", "MLL"]
author_list = ["Eduardo G. Gusmao", "Giuseppe Petrosino", "Argyris Papantonis", "Vassilis Roukos"]
corresponding_mail = "eduardogade@gmail.com"
license_type = "GPL"
data_config_path_file_name = os.path.join(os.getcwd(), "data")
package_data_dictionary = {"src": [path.basename(data_config_path_file_name)]}

# External scripts
external_scripts = []
for tool_option in options.param_tools:
    for e in tools_dictionary[tool_option][3]:
        external_scripts.append(e)

# Fetching Additional Structural Files
readme_file_name = path.join(path.dirname(path.abspath(__file__)), "README.rst")

# Fetching Long Description
readme_file = open(readme_file_name, "r")
long_description = readme_file.read()
readme_file.close()

# Setup Function

setup(name="GotheEtAl",
      version=current_version,
      description=short_description,
      long_description=long_description,
      classifiers=classifiers_list,
      keywords=", ".join(keywords_list),
      author=", ".join(author_list),
      author_email=corresponding_mail,
      license=license_type,
      packages=find_packages(),
      package_data=package_data_dictionary,
      entry_points=current_entry_points,
      install_requires=current_install_requires,
      scripts=external_scripts,
      platforms=supported_platforms)

