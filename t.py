#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import with_statement
from collections import defaultdict
from math import sqrt
from random import shuffle, randint

# 欧几里德距离相关度
def sim_distance(prefs, v1, v2):
    si = {}
    # 找到共同的元素
    for item in prefs[v1]:
        if item in prefs[v2]:
            si[item] = 1

    # 没有共同的元素
    if len(si) == 0: return 0

    # 计算所有元素的差值平方和
    sum_of_squares = sum([pow(prefs[v1][item] - prefs[v2][item], 2) for item in si])

    # 欧几里德相关度
    # 加1是为了避免除数为0的情况
    # 0相关度小，1相关度大
    return 1 / (1 + sqrt(sum_of_squares))

# Pearson相关度
def sim_pearson(prefs, v1, v2):
    si = {}
    # 找到共同的元素
    for item in prefs[v1]:
        if item in prefs[v2]:
            si[item] = 1

    # 没有共同的元素
    n = len(si)
    if n == 0: return 0

    # 计算所有元素的和
    sum1 = sum([prefs[v1][item] for item in si])
    sum2 = sum([prefs[v2][item] for item in si])

    # 计算所有元素的平方和
    sum1Sq = sum([pow(prefs[v1][item], 2) for item in si])
    sum2Sq = sum([pow(prefs[v2][item], 2) for item in si])

    # 计算乘机和
    pSum = sum([prefs[v1][item] * prefs[v2][item] for item in si])

    # Pearson相关度
    num = pSum - (sum1 * sum2 / n)
    den = sqrt((sum1Sq - pow(sum1, 2) / n) * (sum2Sq - pow(sum2, 2) / n))
    if den == 0: return 0

    # -1相关度小，1相关度大
    r = num / den
    return r

# Cosine相关度
def sim_cosine(prefs, v1, v2):
    si = {}
    # 找到共同的元素
    for item in prefs[v1]:
        if item in prefs[v2]:
            si[item] = 1

    # 没有共同的元素
    n = len(si)
    if n == 0: return 0

    # 计算乘机和
    pSum = sum([prefs[v1][item] * prefs[v2][item] for item in si])

    # 计算所有元素的平方和
    sum1Sq = sum([pow(prefs[v1][item], 2) for item in si])
    sum2Sq = sum([pow(prefs[v2][item], 2) for item in si])

    return pSum / (sqrt(sum1Sq) * sqrt(sum2Sq))

# Jaccard Index相关度
def sim_jaccard(prefs, v1, v2):
    si = {}
    # 找到共同的元素
    for item in prefs[v1]:
        if item in prefs[v2]:
            si[item] = 1

    # 没有共同的元素
    n = len(si)
    if n == 0: return 0

    # 计算乘机和
    pSum = sum([prefs[v1][item] * prefs[v2][item] for item in si])

    # 计算所有元素的平方和
    sum1Sq = sum([pow(prefs[v1][item], 2) for item in si])
    sum2Sq = sum([pow(prefs[v2][item], 2) for item in si])

    return pSum / (sum1Sq + sum2Sq - pSum)

# Tanimoto系数
def sim_tanimoto(prefs, v1, v2):

    # 出现在Person1中
    c1 = len(prefs[v1])
    # 出现在Person2中
    c2 = len(prefs[v2])
    # 两个都出现的
    c = [item for item in prefs[v1] if item in prefs[v2]]
    shr = len(c)

    return float(shr) / (c1 + c2 - shr)

# 曼哈顿距离相关度
def sim_taxicab(prefs, v1, v2):
    si = {}
    # 找到共同的元素
    for item in prefs[v1]:
        if item in prefs[v2]:
            si[item] = 1

    # 没有共同的元素
    n = len(si)
    if n == 0: return 0

    # 计算所有元素的差值绝对值和
    sum_of_abs = sum([abs(prefs[v1][item] - prefs[v2][item]) for item in si])

    # 曼哈顿相关度
    # 加1是为了避免除数为0的情况
    # 0相关度小，1相关度大
    return 1 / (1 + sum_of_abs)

# 转换
def transformPrefs(prefs):
    result = defaultdict(dict)
    for person in prefs:
        for item in prefs[person]:
            result[item][person] = prefs[person][item]
    return result

# 相关度Top
def topMatches(prefs, v, n = 5, similarity = sim_pearson):
    # 计算与其他的相关度
    scores = [(similarity(prefs, v, other), other) for other in prefs if other != v]

    scores.sort()
    scores.reverse()
    return scores[0:n]

# 计算每个物品的相关度
def calculateSimilarItems(prefs, n = 10, similarity = sim_pearson):
    results = {}
    # 以物品为中心的偏好矩阵
    itemPrefs = transformPrefs(prefs)
    c = 0
    for item in itemPrefs:
        c += 1
        if c % 100 == 0:
            print '%d / %d' % (c, len(itemPrefs))
        # 找到相似的物品
        scores = topMatches(itemPrefs, item, n = n, similarity = similarity)
        results[item] = scores
    return results

# 获取某个用户的推荐物品
def getRecommendedItems(prefs, itemMatch, user):
    userRatings = prefs[user]
    scores = {}
    totalSim = {}

    # 循环遍历由当前用户评分的物品
    for (item, rating) in userRatings.items():
        # 循环遍历与当前物品相近的物品
        for (similarity, item2) in itemMatch[item]:
            # 如果该用户已经对当前物品做过评价，则将其忽略
            if item2 in userRatings: continue

            # 评价值与相似度的加权之和
            scores.setdefault(item2, 0)
            scores[item2] += similarity * rating

            # 全部相似度之和
            totalSim.setdefault(item2, 0)
            totalSim[item2] += similarity

    # 将每个合计值除以加权和，求出平均值
    ranking = [(score / totalSim[item], item) for item, score in scores.items()]

    ranking.sort()
    ranking.reverse()
    return ranking

def main(path):
    datas = defaultdict(dict)
    with open(path, 'r') as fp:
        for line in fp:
            uid, tid, score = line.strip().split('\t')
            datas[uid][tid] = float(score)
    calculateSimilarItems(datas)
    '''
    for uid in datas:
        for auid in datas:
            if uid != auid:
                simi = sim_distance(datas, uid, auid)
                if (simi != 0):
                    print 'distance', simi, uid, auid
                simi = sim_pearson(datas, uid, auid)
                if (simi != 0):
                    print 'pearson', simi, uid, auid
                simi = sim_cosine(datas, uid, auid)
                if (simi != 0):
                    print 'cosine', simi, uid, auid
                simi = sim_jaccard(datas, uid, auid)
                if (simi != 0):
                    print 'jaccard', simi, uid, auid
                simi = sim_tanimoto(datas, uid, auid)
                if (simi != 0):
                    print 'tanimoto', simi, uid, auid
                simi = sim_taxicab(datas, uid, auid)
                if (simi != 0):
                    print 'taxicab', simi, uid, auid
                    print '================'
    '''

def shuffle_data(path):
    uids = []
    tids = []
    with open(path, 'r') as fp:
        for line in fp:
            uid, tid, score = line.strip().split('\t')
            uids.append(uid)
            tids.append(tid)
    shuffle(uids)
    shuffle(tids)
    for i in range(len(uids)):
        print '%s\t%s\t%s' % (uids[i], tids[i], randint(1, 5))

if __name__ == '__main__':
    #shuffle_data('recommendation.txt')
    main('recommendations.txt')
