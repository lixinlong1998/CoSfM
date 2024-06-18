import numpy as np
import matplotlib.pyplot as plt
from PIL import Image


def create_chessboard_images(img1_path, img2_path, output1_path, output2_path):
    # 读取输入图片
    img1 = np.array(Image.open(img1_path))
    img2 = np.array(Image.open(img2_path))

    # 确保两张图片大小一致
    assert img1.shape == img2.shape, "两张图片大小不一致"

    # 获取图片大小和棋盘尺寸
    height, width, channels = img1.shape
    rows, cols = 8, 8
    block_height = height // rows
    block_width = width // cols

    # 创建空白的棋盘图
    chess1 = np.zeros_like(img1)
    chess2 = np.zeros_like(img2)

    # 填充棋盘图
    for i in range(rows):
        for j in range(cols):
            if (i + j) % 2 == 0:
                chess1[i * block_height:(i + 1) * block_height, j * block_width:(j + 1) * block_width, :] = img1[
                                                                                                            i * block_height:(
                                                                                                                                     i + 1) * block_height,
                                                                                                            j * block_width:(
                                                                                                                                    j + 1) * block_width,
                                                                                                            :]
                chess2[i * block_height:(i + 1) * block_height, j * block_width:(j + 1) * block_width, :] = img2[
                                                                                                            i * block_height:(
                                                                                                                                     i + 1) * block_height,
                                                                                                            j * block_width:(
                                                                                                                                    j + 1) * block_width,
                                                                                                            :]
            else:
                chess1[i * block_height:(i + 1) * block_height, j * block_width:(j + 1) * block_width, :] = img2[
                                                                                                            i * block_height:(
                                                                                                                                     i + 1) * block_height,
                                                                                                            j * block_width:(
                                                                                                                                    j + 1) * block_width,
                                                                                                            :]
                chess2[i * block_height:(i + 1) * block_height, j * block_width:(j + 1) * block_width, :] = img1[
                                                                                                            i * block_height:(
                                                                                                                                     i + 1) * block_height,
                                                                                                            j * block_width:(
                                                                                                                                    j + 1) * block_width,
                                                                                                            :]

    # 保存输出图片
    Image.fromarray(chess1).save(output1_path)
    Image.fromarray(chess2).save(output2_path)
    print('finished.')


# 示例用法
img1_path = r'D:\Research\20221223_CoSfM\Figure\Fig11_GroundTruth\_Be1e3_unstable2_dom1.jpg'
img2_path = r'D:\Research\20221223_CoSfM\Figure\Fig11_GroundTruth\_Be1e3_unstable2_dom2.jpg'
output1_path = img1_path.replace('_dom1.jpg', '_chess1.jpg')
output2_path = img2_path.replace('_dom2.jpg', '_chess2.jpg')
create_chessboard_images(img1_path, img2_path, output1_path, output2_path)
