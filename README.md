# CoSfM: CFTM
A tool for high precision co-registration of multi-epoch images photogrammetry achieved with Agisoft Metashape Python API and COLMAP (optional).

# CFTM User Manual Version 1.0
update 2024.06.24

# About

CFTM is a type of co-SfM method for co-registering multi-epoch aerial images while achieving photogrammetry simultaneously. With the novel matching strategy - 'matching on feature tracks' and designed optimization algorithms, it can generate dense and evenly distributed common tie points (CTPs), which may not insufficient in traditional incremental SfM, therefore achieving high-precision co-alignment accuracy. Apart from that, the part iterative optimization of CFTM can improve the accuracy of traditional co-alignment in normal cases. The output of CFTM is the refined image poses, witch means the later generated products (such as DEMs, DOMs, and 3D meshs) are directly located in the same coordinate system with nearly unbias, enabling high-precision spatiotemporal analysis (such as surface deformation monitoring).
![image](https://github.com/lixinlong1998/CoSfM/assets/73974003/50e69665-5896-4140-9630-a3a19e36b490)

# Download
Source code is available at [https://github.com/lixinlong1998/CoSfM](https://github.com/lixinlong1998/CoSfM).

Documentation (English version) is available at [https://blog.csdn.net/LXLng/article/details/134613468](https://blog.csdn.net/LXLng/article/details/134613468).

Documentation (Chinese version) is available at [https://blog.csdn.net/LXLng/article/details/136598055](https://blog.csdn.net/LXLng/article/details/136598055).

Example dataset and Tutorial is available from BaiduNetDisk: [https://pan.baidu.com/s/1pQaQGsnx3l-Th9ohO8zUGQ?pwd=yieb](https://pan.baidu.com/s/1pQaQGsnx3l-Th9ohO8zUGQ?pwd=yieb)

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

## Step 4  Install COLMAP (Recommended)

> This script employ COLMAP version 3.6, and other versions have not been tested.

Official installation tutorial: [Official Website](https://demuc.de/colmap/), [GitHub Releases](https://github.com/colmap/colmap/releases)

Other installation tutorials: [Installation, Debugging, 3D Reconstruction Practice, and Intermediate Results Output of Colmap under Windows 10](https://blog.csdn.net/LXLng/article/details/121039782)

Download COLMAP 3.6 release at [https://github.com/colmap/colmap/releases/tag/3.6](https://github.com/colmap/colmap/releases/tag/3.6)

![image](https://github.com/lixinlong1998/CoSfM/assets/73974003/f32e5459-95f3-4903-bc8d-caab0d353fee)

Click **COLMAP-3.6-windows-cuda.zip** or **COLMAP-3.6-windows-no-cuda.zip** (depending on whether your device has NVIDIA GPU enabled) to download. After downloading, you will get a compressed file. Simply unzip it. Next, double-click "RUN_TESTS" to test the environment. If there are no issues, the COLMAP are successfully installed and you can directly double-click "COLMAP" to enter the software's graphical interface. Please remember the path to COLMAP.bat, as it will be needed for subsequent configurations.

#  Tutorial

After installation, please check the documentation to get the configuration and parameter settings

#  License
Please note the licenses of COLMAP if you choose to use it


