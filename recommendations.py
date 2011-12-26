#!/usr/bin/env python
# -*- coding: utf-8 -*-

from math import sqrt

critics={'Lisa Rose': {'Lady in the Water': 2.5, 'Snakes on a Plane': 3.5,
    'Just My Luck': 3.0, 'Superman Returns': 3.5, 'You, Me and Dupree': 2.5,
    'The Night Listener': 3.0},
    'Gene Seymour': {'Lady in the Water': 3.0, 'Snakes on a Plane': 3.5,
        'Just My Luck': 1.5, 'Superman Returns': 5.0, 'The Night Listener': 3.0,
        'You, Me and Dupree': 3.5},
    'Michael Phillips': {'Lady in the Water': 2.5, 'Snakes on a Plane': 3.0,
        'Superman Returns': 3.5, 'The Night Listener': 4.0},
    'Claudia Puig': {'Snakes on a Plane': 3.5, 'Just My Luck': 3.0,
        'The Night Listener': 4.5, 'Superman Returns': 4.0,
        'You, Me and Dupree': 2.5},
    'Mick LaSalle': {'Lady in the Water': 3.0, 'Snakes on a Plane': 4.0,
        'Just My Luck': 2.0, 'Superman Returns': 3.0, 'The Night Listener': 3.0,
        'You, Me and Dupree': 2.0},
    'Jack Matthews': {'Lady in the Water': 3.0, 'Snakes on a Plane': 4.0,
        'The Night Listener': 3.0, 'Superman Returns': 5.0, 'You, Me and Dupree': 3.5},
    'Toby': {'Snakes on a Plane':4.5,'You, Me and Dupree':1.0,'Superman Returns':4.0}}

# 欧几里德距离相关度
def sim_distance(prefs, person1, person2):
    si = {}
    # 找到共同的元素
    for item in prefs[person1]:
        if item in prefs[person2]:
            si[item] = 1

    # 没有共同的元素
    if len(si) == 0: return 0

    # 计算所有元素的差值平方和
    # sum_of_squares = sum([pow(prefs[person1][item] - prefs[person2][item], 2) for item in prefs[person1] if item in prefs[person2]])
    sum_of_squares = sum([pow(prefs[person1][item] - prefs[person2][item], 2) for item in si])

    # 欧几里德相关度
    # 加1是为了避免除数为0的情况
    # 0相关度小，1相关度大
    return 1 / (1 + sqrt(sum_of_squares))

# Pearson相关度
def sim_pearson(prefs, person1, person2):
    si = {}
    # 找到共同的元素
    for item in prefs[person1]:
        if item in prefs[person2]:
            si[item] = 1

    # 没有共同的元素
    n = len(si)
    if n == 0: return 0

    # 计算所有元素的和
    sum1 = sum([prefs[person1][item] for item in si])
    sum2 = sum([prefs[person2][item] for item in si])

    # 计算所有元素的平方和
    sum1Sq = sum([pow(prefs[person1][item], 2) for item in si])
    sum2Sq = sum([pow(prefs[person2][item], 2) for item in si])

    # 计算乘机和
    pSum = sum([prefs[person1][item] * prefs[person2][item] for item in si])

    # Pearson相关度
    num = pSum - (sum1 * sum2 / n)
    den = sqrt((sum1Sq - pow(sum1, 2) / n) * (sum2Sq - pow(sum2, 2) / n))
    if den == 0: return 0

    # -1相关度小，1相关度大
    r = num / den
    return r

# Cosine相关度
def sim_cosine(prefs, person1, person2):
    si = {}
    # 找到共同的元素
    for item in prefs[person1]:
        if item in prefs[person2]:
            si[item] = 1

    # 没有共同的元素
    n = len(si)
    if n == 0: return 0

    # 计算乘机和
    pSum = sum([prefs[person1][item] * prefs[person2][item] for item in si])

    # 计算所有元素的平方和
    sum1Sq = sum([pow(prefs[person1][item], 2) for item in si])
    sum2Sq = sum([pow(prefs[person2][item], 2) for item in si])

    return pSum / (sqrt(sum1Sq) * sqrt(sum2Sq))

# Jaccard Index相关度
def sim_jaccard(prefs, person1, person2):
    si = {}
    # 找到共同的元素
    for item in prefs[person1]:
        if item in prefs[person2]:
            si[item] = 1

    # 没有共同的元素
    n = len(si)
    if n == 0: return 0

    # 计算乘机和
    pSum = sum([prefs[person1][item] * prefs[person2][item] for item in si])

    # 计算所有元素的平方和
    sum1Sq = sum([pow(prefs[person1][item], 2) for item in si])
    sum2Sq = sum([pow(prefs[person2][item], 2) for item in si])

    return pSum / (sum1Sq + sum2Sq - pSum)

# Tanimoto系数
def sim_tanimoto(prefs, person1, person2):

    # 出现在Person1中
    c1 = len(prefs[person1])
    # 出现在Person2中
    c2 = len(prefs[person2])
    # 两个都出现的
    c = [item for item in prefs[person1] if item in prefs[person2]]
    shr = len(c)

    return float(shr) / (c1 + c2 - shr)

# 曼哈顿距离相关度
def sim_taxicab(prefs, person1, person2):
    si = {}
    # 找到共同的元素
    for item in prefs[person1]:
        if item in prefs[person2]:
            si[item] = 1

    # 没有共同的元素
    n = len(si)
    if n == 0: return 0

    # 计算所有元素的差值绝对值和
    # sum_of_abs = sum([abs(prefs[person1][item] - prefs[person2][item]) for item in prefs[person1] if item in prefs[person2]])
    sum_of_abs = sum([abs(prefs[person1][item] - prefs[person2][item]) for item in si])

    # 曼哈顿相关度
    # 加1是为了避免除数为0的情况
    # 0相关度小，1相关度大
    return 1 / (1 + sum_of_abs)

if __name__ == '__main__':
    print(sim_distance(critics, 'Lisa Rose', 'Mick LaSalle'))
    print(sim_distance(critics, 'Toby', 'Claudia Puig'))
    print(sim_pearson(critics, 'Lisa Rose', 'Mick LaSalle'))
    print(sim_pearson(critics, 'Toby', 'Claudia Puig'))
    print(sim_cosine(critics, 'Lisa Rose', 'Mick LaSalle'))
    print(sim_cosine(critics, 'Toby', 'Claudia Puig'))
    print(sim_jaccard(critics, 'Lisa Rose', 'Mick LaSalle'))
    print(sim_jaccard(critics, 'Toby', 'Claudia Puig'))
    print(sim_tanimoto(critics, 'Lisa Rose', 'Mick LaSalle'))
    print(sim_tanimoto(critics, 'Toby', 'Claudia Puig'))
    print(sim_taxicab(critics, 'Lisa Rose', 'Mick LaSalle'))
    print(sim_taxicab(critics, 'Toby', 'Claudia Puig'))
