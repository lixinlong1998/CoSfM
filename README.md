# CoSfM: CFTM
A tool for high precision co-registration of multi-epoch images photogrammetry based on Agisoft Metashape Python API and COLMAP.

# CFTM User Manual Version 1.0
update 2024.03.010

# About

CFTM is a tool developed based on Agisoft Metashape. It is inspired by advanced co-alignment algorithms and improves the matching method and adds filtering optimization algorithms. It can achieve high-precision registration and photogrammetric processing of multi-period aerial images simultaneously, placing the generated multi-period three-dimensional spatial data directly in the same reference coordinate system, enabling high-precision spatiotemporal analysis applications, such as surface deformation monitoring.

# Download
The source code is available at [https://github.com/lixinlong1998/CoSfM](https://github.com/lixinlong1998/CoSfM).

The English version of documentation is available at [https://blog.csdn.net/LXLng/article/details/134613468](https://blog.csdn.net/LXLng/article/details/134613468).

The Chinese version of documentation is available at [https://blog.csdn.net/LXLng/article/details/136598055](https://blog.csdn.net/LXLng/article/details/136598055).

The dataset for test is available at [https://pan.baidu.com/s/1nNklV2M1Zs02UkWtGG_tLw](https://pan.baidu.com/s/1nNklV2M1Zs02UkWtGG_tLw)

Alternatively, you can download all properties mentioned above from BaiduNetDisk link:
[https://pan.baidu.com/s/1pQaQGsnx3l-Th9ohO8zUGQ?pwd=yieb](https://pan.baidu.com/s/1pQaQGsnx3l-Th9ohO8zUGQ?pwd=yieb)

# Installation

CFTM v1.0 is developed using the Python API of Agisoft Metashape Pro v1.8.5 in a Windows environment. Since Metashape Pro v1.8.5 does not provide interfaces for feature point extraction and raw feature matching, the **Feature Point Extraction Module** utilizes the feature extraction function (SIFT-GPU) from the third-party open-source library COLMAP, and the **Result Analysis Module** utilizes the feature matching function from COLMAP. In addition, some functionalities depend on third-party library functions. Although Metashape Pro v1.8.5 comes with a pre-installed Python 3.8 environment with some third-party libraries, additional installation and/or upgrading of some third-party libraries are still required.


## Step 1  Installation of Agisoft Metashape Pro v1.8.5

> If you are new to Metashape, please note that the software is paid and requires a license purchase. Agisoft offers an Educational License version, details of which can be found [here](https://www.agisoft.com/buy/licensing-options/).
> Below are some related links to help you quickly understand and reference Metashape:
> [Official Website](https://www.agisoft.com/downloads/installer/), [Beginner's Tutorial](https://www.agisoft.com/support/tutorials/), [Forum](https://www.agisoft.com/forum/), [Knowledge Base](https://agisoft.freshdesk.com/support/solutions)

* [ ] Visit the [Agisoft official website](https://www.agisoft.com/downloads/installer/) and download the installer for Metashape Pro v1.8.5. Note that the official website may remove previous versions with version updates. It is recommended for users to search for historical versions through Google (for everyone's convenience, I will post the download link for the installation package in the comments).
* [ ] Run the downloaded installer and follow the installation wizard instructions to complete the installation process.
* [ ] During the installation process, you may need to accept the license agreement, choose the installation path, and select other options.
## Step 2  Adding Metashape.exe to the Environment Variables

* [ ] Right-click on `This PC` and select `Properties`, then find `Advanced system settings` and open the panel.

![1](https://img-blog.csdnimg.cn/direct/a5d8910f703743b88dcf12cc68048c8d.png#pic_center)
* [ ] Open the Environment Variables panel, select `Path` under `System variables`, and click `Edit`.

![2](https://img-blog.csdnimg.cn/direct/cff55c0155424fa78c429a11842b9170.png#pic_center)

* [ ] Click `New`, then enter the installation path (the folder path containing Metashape.exe), and proceed to `OK` all panels to complete the environment variable setup.

![3](https://img-blog.csdnimg.cn/direct/65e75d1da67b4fc5af28d66707fb463b.png#pic_center)

![4](https://img-blog.csdnimg.cn/direct/afd77a456d5147e0b5f291ca37e12211.png#pic_center)

## Step 3  Installing Third-Party Python Libraries

> Metashape 1.8.5 comes with a Python 3.8 environment and some basic Python libraries. Our plugin is developed under the Python 3.8 environment provided by Metashape 1.8.5. However, since some complex functionalities require additional library functions to be implemented, it is necessary to install new third-party library functions in the existing Python 3.8 environment.

The official tutorial is available here: [How to install external Python module to Metashape Professional package](https://agisoft.freshdesk.com/support/solutions/articles/31000136860-how-to-install-external-python-module-to-metashape-professional-package). We provide 2 methods to install third-party library functions in the specified Python environment on Windows, using numpy as an example:

**Method 1: Automatic Download**

* [ ] Press Win+R, type `cmd`, and press Enter to open the command line cmd.exe.
* [ ] Type `cd /d <G:\UAV_SOFTWARE\Metashape\Metashape1.8.5\python>` (replace with your installation path).
* [ ] Type `python -m pip install numpy` to install the numpy library.
* [ ] If the library is not found, you can try specifying the source URL afterward. For example, to use the Tsinghua source, type `python -m pip install numpy -i https://pypi.tuna.tsinghua.edu.cn/simple/`.

**Method 2: Manual Download**

* [ ] Download the appropriate version of the wheel file from the source website in advance, e.g., `numpy-1.24.2-cp38-cp38-win_amd64.whl`.
* [ ] Press Win+R, type `cmd`, and press Enter to open the command line cmd.exe.
* [ ] Type `cd /d <G:\UAV_SOFTWARE\Metashape\Metashape1.8.5\python>`, replacing with your installation path.
* [ ] Execute the command `python -m pip install <wheel_path>`.

Below are the third-party library functions that need to be installed:

> 1. The versions of these libraries have been tested to be suitable for Metashape Pro v1.8.5.
> 2. It is recommended to manually download the wheel file for GDAL and then use Method 2 for installation.

|       Library Name      |              Version or Wheel               | Recommended Method |
| :------------------------: | :------------------------------------------: | :---------------------: |
|            numpy            |                    1.23.5                    |         Automatic Download         |
|            scipy            |                    1.10.1                    |         Automatic Download         |
|           pandas           |                     2.0.0                     |         Automatic Download         |
|          networkx          |                       3.1                       |         Automatic Download         |
|            GDAL            | GDAL‑3.4.3‑cp38‑cp38‑win_amd64.whl |              [Manual Download](https://www.lfd.uci.edu/~gohlke/pythonlibs/#pygame)               |
|        opencv-python       |                     4.7                      |         Automatic Download         |
| opencv-contrib-python |                     4.7                      |         Automatic Download         |
|         matplotlib         |                    3.7.1                     |         Automatic Download         |
|       python-dateutil      |                    2.8.2                     |         Automatic Download         |
|            pytz            |                    2023.3                    |         Automatic Download         |

Installation Process:

* [ ] First, download the GDAL installation package, and then execute the following commands line by line.
* [ ] If you encounter a corrupted pip issues during the installation process, you can follow this method  [^2].
[^2]:Regarding the repair method for a broken pip in the Python environment bundled with Metashape:
Run the following commands line by line in the command line, **making sure to replace the paths with your own**:`cd /d <G:\UAV_SOFTWARE\Metashape\Metashape1.8.5\python>`, `python <...CFTM_v1.0\toolbox\get-pip.py>`, `python -m pip install --upgrade pip`, After executing these commands, pip should function properly. However, upon reopening Metashape, it will automatically update and restore the original pip. If you wish to prevent this, you can try turning off the automatic update option in Metashape (*Tools → Preferences → Miscellaneous → Check for updates on program start*).


```bash
cd /d <G:\UAV_SOFTWARE\Metashape\Metashape1.8.5\python>
python -m pip install --upgrade pip
python -m pip install numpy==1.23.5
python -m pip install scipy==1.10.1
python -m pip install pandas==2.0.0
python -m pip install matplotlib==3.7.1
python -m pip install networkx==3.1
python -m pip install python-dateutil==2.8.2
python -m pip install pytz==2023.3
python -m pip install opencv-python==4.7.0.72 -i https://pypi.tuna.tsinghua.edu.cn/simple/
python -m pip install opencv-contrib-python==4.7.0.72 -i https://pypi.tuna.tsinghua.edu.cn/simple/
python -m pip install G:\UAV_SOFTWARE\Metashape\Python3Module\GDAL-3.4.3-cp38-cp38-win_amd64.whl
```

* [ ] After installation, type the following commands line by line in the command prompt to check the libraries installed in the current Metashape Python environment.

```bash
cd /d <G:\UAV_SOFTWARE\Metashape\Metashape1.8.5\python>
python -m pip list
```


| Library | Version | Library | Version |
| ----------------------- | ---------- | ---------------------- | ---------- |
| attrs                 | 23.1.0   | backcall             | 0.2.0    |
| click                 | 8.1.3    | colorama             | 0.4.3    |
| ConfigArgParse        | 1.5.3    | contourpy            | 1.0.7    |
| cycler                | 0.11.0   | dash                 | 2.9.3    |
| dash-core-components  | 2.0.0    | dash-html-components | 2.0.0    |
| dash-table            | 5.0.0    | decorator            | 4.4.2    |
| fastjsonschema        | 2.16.3   | Flask                | 2.2.3    |
| fonttools             | 4.39.3   | GDAL                 | 3.4.3    |
| importlib-metadata    | 6.5.0    | importlib-resources  | 5.12.0   |
| ipykernel             | 5.3.4    | ipython              | 7.16.1   |
| ipython-genutils      | 0.2.0    | ipywidgets           | 8.0.6    |
| itsdangerous          | 2.1.2    | jedi                 | 0.17.2   |
| Jinja2                | 3.1.2    | jsonschema           | 4.17.3   |
| jupyter-client        | 6.1.6    | jupyter-core         | 4.6.3    |
| jupyterlab-widgets    | 3.0.7    | kiwisolver           | 1.4.4    |
| llvmlite              | 0.39.1   | MarkupSafe           | 2.1.2    |
| matplotlib            | 3.7.1    | nbformat             | 5.7.0    |
| networkx              | 3.1      | numba                | 0.56.4   |
| numpy                 | 1.23.5   | open3d               | 0.17.0   |
| opencv-contrib-python | 4.7.0.72 | opencv-python        | 4.7.0.72 |
| packaging             | 23.1     | pandas               | 2.0.0    |
| parso                 | 0.7.1    | pexpect              | 4.8.0    |
| pickleshare           | 0.7.5    | Pillow               | 9.5.0    |
| pip                   | 24.0     | pkgutil_resolve_name | 1.3.10   |
| plotly                | 5.14.1   | prompt-toolkit       | 3.0.5    |
| ptyprocess            | 0.6.0    | Pygments             | 2.6.1    |
| pyparsing             | 3.0.9    | pyrsistent           | 0.19.3   |
| PySide2               | 5.15.2.1 | python-dateutil      | 2.8.1    |
| pytz                  | 2023.3   | pywin32              | 228      |
| pyzmq                 | 19.0.1   | qtconsole            | 4.7.5    |
| QtPy                  | 1.9.0    | scipy                | 1.10.1   |
| setuptools            | 49.2.0   | shiboken2            | 5.15.2.1 |
| shiboken2-generator   | 5.15.2.1 | six                  | 1.15.0   |
| tenacity              | 8.2.2    | tornado              | 6.0.4    |
| traitlets             | 4.3.3    | tzdata               | 2023.3   |
| wcwidth               | 0.2.5    | Werkzeug             | 2.2.3    |
| wheel                 | 0.29.0   | widgetsnbextension   | 4.0.7    |
| zipp                  | 3.15.0   |                      |          |


## Step 4  Install COLMAP (Recommended)

> This script employ COLMAP version 3.6, and other versions have not been tested.

Official installation tutorial: [Official Website](https://demuc.de/colmap/), [GitHub Releases](https://github.com/colmap/colmap/releases)
Other installation tutorials: [Installation, Debugging, 3D Reconstruction Practice, and Intermediate Results Output of Colmap under Windows 10](https://blog.csdn.net/LXLng/article/details/121039782)

Download COLMAP 3.6 release at [https://github.com/colmap/colmap/releases/tag/3.6](https://github.com/colmap/colmap/releases/tag/3.6)
![image](https://github.com/lixinlong1998/CoSfM/assets/73974003/f32e5459-95f3-4903-bc8d-caab0d353fee)
Click **COLMAP-3.6-windows-cuda.zip** or **COLMAP-3.6-windows-no-cuda.zip** (Depending on whether your device has NVIDIA GPU enabled.) to download
After downloading, you will get a compressed file. Simply unzip it.
"bin" contains the executable binary files of the software.
"lib" contains third-party library files used.
"COLMAP" is the program. Simply double-click to run.
"RUN_TESTS" is the testing file.
Next, double-click "RUN_TESTS" to test the environment. If there are no issues, the COLMAP are successfully installed and you can directly double-click "COLMAP" to enter the software's graphical interface. Please remember the path to COLMAP.bat, as it will be needed for subsequent configurations.

#  License


