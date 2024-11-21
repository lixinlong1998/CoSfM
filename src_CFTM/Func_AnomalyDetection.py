'''
Q:如何使用高斯混合模型（Gaussian Mixture Model，GMM）分解混合高斯分布？
ChatGPT:
初始化模型参数：随机初始化两个高斯分布的均值和方差，以及混合系数。
    Expectation步骤：对于每个数据点，计算其属于每个高斯分布的后验概率，即它来自第一个高斯分布的概率和它来自第二个高斯分布的概率。
    Maximization步骤：根据每个数据点的后验概率重新计算高斯分布的参数，即更新均值和方差，以及混合系数。
重复步骤2和3，直到收敛。
收敛后，我们就可以得到两个高斯分布的均值和方差，以及它们的混合系数。由于我们已知两个高斯分布的均值和方差的大小关系，可以根据这个信息对它们进行标识。
注意，GMM是一种迭代算法，初始参数的选择可能会影响最终结果。因此，可能需要尝试多组不同的初始参数，并选择最终结果最好的一组。

如果您不想使用scikit-learn库，可以使用scipy.stats中的norm函数来拟合高斯分布，并使用EM算法实现GMM。具体来说，可以使用scipy.optimize中的
minimize函数来最大化似然函数，并迭代更新参数。
'''
import numpy as np
from scipy.stats import norm
from scipy.optimize import minimize


def GMM(data, means=(0, 0), stds=(1, 1), weights=(0.5, 0.5)):
    '''
    input:
        data = np.array([item,...])
        means=(mean1,mean2,...)
        stds=(std1,std2,...)
        weights=(w1,w2,...), where sum(weights)==1
    output:
        data = list[np.array,np.array,...]
        params = np.array([[mean,std,weight],[mean,std,weight],...])
    '''

    def log_likelihood(params, data):
        # 定义似然函数
        n = len(data)
        k = len(means)
        means = params[:k]
        stds = params[k:2 * k]
        weights = params[2 * k:]
        ll = 0
        for i in range(n):
            p = 0
            for j in range(k):
                p += weights[j] * norm.pdf(data[i], loc=means[j], scale=stds[j])
            ll += np.log(p)
        return -ll

    # 初始化参数
    params = np.concatenate([means, stds, weights])

    # 最大化似然函数
    result = minimize(log_likelihood, params, args=(data,))
    means = result.x[:2]
    stds = result.x[2:4]
    weights = result.x[4:]

    # 确定高斯分布的顺序
    if means[0] > means[1]:
        means = np.flip(means)
        stds = np.flip(stds)
        weights = np.flip(weights)
    print("高斯分布1：均值={:.2f}，方差={:.2f}，权重={:.2f}".format(means[0], stds[0] ** 2, weights[0]))
    print("高斯分布2：均值={:.2f}，方差={:.2f}，权重={:.2f}".format(means[1], stds[1] ** 2, weights[1]))



'''
Q:如何判断列表中每个数据属于哪个高斯分布?
ChatGPT:在使用高斯混合模型对数据进行建模后，我们可以通过计算每个数据点对应于每个高斯分布的后验概率来确定该数据点属于哪个分布。具体来说，我们可以使用
scipy.stats.norm中的pdf函数计算每个数据点属于每个高斯分布的概率密度函数值，然后根据每个高斯分布的权重将其归一化得到后验概率。
'''
import numpy as np
from scipy.stats import norm

# 假设我们已经拟合了高斯混合模型，并得到了每个高斯分布的均值、方差和权重
means = [1, 10]
stds = [1, 5]
weights = [0.4, 0.6]

# 假设我们有一个包含一些数据点的列表
data = np.array([2, 3, 4, 11, 12, 13])

# 计算每个数据点属于每个高斯分布的后验概率
n = len(data)
k = len(means)
posterior_probs = np.zeros((n, k))
for i in range(n):
    for j in range(k):
        posterior_probs[i, j] = weights[j] * norm.pdf(data[i], loc=means[j], scale=stds[j])
    posterior_probs[i] /= np.sum(posterior_probs[i])

# 打印结果
for i in range(n):
    print("数据点{}属于高斯分布1的概率为{:.2f}，属于高斯分布2的概率为{:.2f}".format(data[i], posterior_probs[i, 0], posterior_probs[i, 1]))
