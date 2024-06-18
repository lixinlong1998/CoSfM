import cv2
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image


def compute_and_plot_optical_flow_0(img1_path, img2_path, save_path, base_image_path):
    # 读取输入图片
    img1 = cv2.imread(img1_path)
    img2 = cv2.imread(img2_path)
    base = cv2.imread(base_image_path)
    # 检查图片大小是否一致
    assert img1.shape == img2.shape, "两张图片大小不一致"

    # 转换为灰度图
    gray1 = cv2.cvtColor(img1, cv2.COLOR_BGR2GRAY)
    gray2 = cv2.cvtColor(img2, cv2.COLOR_BGR2GRAY)
    base_grey = cv2.cvtColor(base, cv2.COLOR_BGR2GRAY)

    # 计算光流
    flow = cv2.calcOpticalFlowFarneback(gray1, gray2, None, 0.5, 10, 50, 3, 7, 1.2, 0)

    # 获取光流的水平和垂直分量
    flow_x = flow[..., 0]
    flow_y = flow[..., 1]

    # 绘制光流场
    interval = 200
    plt.figure(figsize=(10, 10), dpi=300)
    plt.imshow(base_grey, cmap='gray')
    plt.axis('off')  # 关闭坐标轴
    h, w = gray1.shape
    y, x = np.mgrid[0:h:interval, 0:w:interval]
    # scale: 箭头的缩放比例，它决定了箭头的长度。在这里，我们将其设置为interval*ratio，ratio越大，箭头越短
    plt.quiver(x, y, flow_x[::interval, ::interval], flow_y[::interval, ::interval], color='red',
               scale=interval * 5,
               angles='xy')
    plt.savefig(save_path, bbox_inches='tight', pad_inches=0)
    # plt.show()


import cv2
import numpy as np
import matplotlib.pyplot as plt


def compute_and_plot_optical_flow2(img1_path, img2_path, save_path):
    # 读取输入图片
    img1 = cv2.imread(img1_path)
    img2 = cv2.imread(img2_path)

    # 检查图片大小是否一致
    assert img1.shape == img2.shape, "两张图片大小不一致"

    # 转换为灰度图
    gray1 = cv2.cvtColor(img1, cv2.COLOR_BGR2GRAY)
    gray2 = cv2.cvtColor(img2, cv2.COLOR_BGR2GRAY)

    # 计算光流
    flow = cv2.calcOpticalFlowFarneback(gray1, gray2, None, 0.5, 10, 50, 3, 7, 1.2, 0)

    # 获取图片的宽高
    h, w = gray1.shape
    vis = img1.copy()

    # 设置参数
    step = 100  # 步长
    arrow_scale = 5  # 箭头的缩放比例
    arrow_thickness = 5  # 箭头的厚度
    tip_length = 0.3  # 箭头的尖端长度

    # 绘制箭头
    for y in range(0, h, step):
        for x in range(0, w, step):
            dx, dy = flow[y, x]
            cv2.arrowedLine(vis, (x, y), (int(x + dx * arrow_scale), int(y + dy * arrow_scale)), (0, 0, 255),
                            arrow_thickness, tipLength=tip_length)

    # 绘制图例
    legend_x = 50
    legend_y = 50
    legend_length = 10  # 图例中表示1个像素位移的长度
    cv2.arrowedLine(vis, (legend_x, legend_y), (legend_x + int(legend_length * arrow_scale), legend_y), (0, 0, 255),
                    arrow_thickness, tipLength=tip_length)
    cv2.putText(vis, '1 pixel displacement', (legend_x, legend_y - 10), cv2.FONT_HERSHEY_SIMPLEX, 0.5, (0, 0, 0), 1)

    # 保存和显示结果
    cv2.imwrite(save_path, vis)
    cv2.imshow('Optical Flow', vis)
    cv2.waitKey(0)
    cv2.destroyAllWindows()


# 示例用法
image_path_list = [(r'D:\Research\20221223_CoSfM\Figure\Fig11_GroundTruth\optical flow\_Be1e2_unstable1_dsm1.jpg',
                    r'D:\Research\20221223_CoSfM\Figure\Fig11_GroundTruth\optical flow\_Be1e2_unstable1_dsm2.jpg'),
                   (r'D:\Research\20221223_CoSfM\Figure\Fig11_GroundTruth\optical flow\_Be1e2_unstable2_dsm1.jpg',
                    r'D:\Research\20221223_CoSfM\Figure\Fig11_GroundTruth\optical flow\_Be1e2_unstable2_dsm2.jpg'),
                   (r'D:\Research\20221223_CoSfM\Figure\Fig11_GroundTruth\optical flow\_Be1e3_unstable1_dsm1.jpg',
                    r'D:\Research\20221223_CoSfM\Figure\Fig11_GroundTruth\optical flow\_Be1e3_unstable1_dsm2.jpg'),
                   (r'D:\Research\20221223_CoSfM\Figure\Fig11_GroundTruth\optical flow\_Be1e3_unstable2_dsm1.jpg',
                    r'D:\Research\20221223_CoSfM\Figure\Fig11_GroundTruth\optical flow\_Be1e3_unstable2_dsm2.jpg'),
                   (r'D:\Research\20221223_CoSfM\Figure\Fig11_GroundTruth\optical flow\_Xe2e3_unstable1_dsm1.jpg',
                    r'D:\Research\20221223_CoSfM\Figure\Fig11_GroundTruth\optical flow\_Xe2e3_unstable1_dsm2.jpg'),
                   (r'D:\Research\20221223_CoSfM\Figure\Fig11_GroundTruth\optical flow\_Xe2e3_unstable2_dsm1.jpg',
                    r'D:\Research\20221223_CoSfM\Figure\Fig11_GroundTruth\optical flow\_Xe2e3_unstable2_dsm2.jpg')]
base_path_list = [r'D:\Research\20221223_CoSfM\Figure\Fig11_GroundTruth\_Be1e2_unstable1_dom1.jpg',
                  r'D:\Research\20221223_CoSfM\Figure\Fig11_GroundTruth\_Be1e2_unstable2_dom2.jpg',
                  r'D:\Research\20221223_CoSfM\Figure\Fig11_GroundTruth\_Be1e3_unstable1_dom1.jpg',
                  r'D:\Research\20221223_CoSfM\Figure\Fig11_GroundTruth\_Be1e3_unstable2_dom1.jpg',
                  r'D:\Research\20221223_CoSfM\Figure\Fig11_GroundTruth\_Xe2e3_unstable1_dom1.jpg',
                  r'D:\Research\20221223_CoSfM\Figure\Fig11_GroundTruth\_Xe2e3_unstable2_dom1.jpg']
#
for i, image_path_pair in enumerate(image_path_list):
    img1_path = image_path_pair[0]
    img2_path = image_path_pair[1]
    base_path = base_path_list[i]
    save_path = img1_path.replace('_dsm1.jpg', '_OpticalFlow_red.jpg')
    # save_path = img1_path.replace('_dom1.jpg', '_OpticalFlow.jpg')

    compute_and_plot_optical_flow_0(img1_path, img2_path, save_path, base_path)
    # compute_and_plot_optical_flow2(img1_path, img2_path, save_path)
