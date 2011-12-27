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

# 相关度Top
def topMatches(prefs, person, n = 5, similarity = sim_pearson):
    # 计算与其他人的相关度
    scores = [(similarity(prefs, person, other), other) for other in prefs if other != person]

    scores.sort()
    scores.reverse()
    return scores[0:n]

# 获取推荐
def getRecommendations(prefs, person, similarity = sim_pearson):
    totals = {}
    simSums = {}
    # 遍历所有人的数据
    for other in prefs:
        # 不包括自己
        if other == person: continue
        # 相似度
        sim = similarity(prefs, person, other)
        # 忽略小于等于0的情况
        if sim <= 0: continue
        for item in prefs[other]:
            # 未评价过的
            if item not in prefs[person] or prefs[person][item] == 0:
                # 相似度 * 评价值
                totals.setdefault(item, 0)
                totals[item] += prefs[other][item] * sim
                # 相似度之和
                simSums.setdefault(item, 0)
                simSums[item] += sim

    # 归一化列表
    ranking = [(total / simSums[item], item) for item, total in totals.items()]
    ranking.sort()
    ranking.reverse()
    return ranking

# 转换
def transformPrefs(prefs):
    result = {}
    for person in prefs:
        for item in prefs[person]:
            result.setdefault(item, {})
            result[item][person] = prefs[person][item]
    return result

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

if __name__ == '__main__':
    print "Lisa Rosa and Mick LaSalle similarity"
    print(sim_distance(critics, 'Lisa Rose', 'Mick LaSalle'))
    print(sim_pearson(critics, 'Lisa Rose', 'Mick LaSalle'))
    print(sim_cosine(critics, 'Lisa Rose', 'Mick LaSalle'))
    print(sim_jaccard(critics, 'Lisa Rose', 'Mick LaSalle'))
    print(sim_tanimoto(critics, 'Lisa Rose', 'Mick LaSalle'))
    print(sim_taxicab(critics, 'Lisa Rose', 'Mick LaSalle'))
    print "Toby and Claudia Puig similarity"
    print(sim_distance(critics, 'Toby', 'Claudia Puig'))
    print(sim_pearson(critics, 'Toby', 'Claudia Puig'))
    print(sim_cosine(critics, 'Toby', 'Claudia Puig'))
    print(sim_jaccard(critics, 'Toby', 'Claudia Puig'))
    print(sim_tanimoto(critics, 'Toby', 'Claudia Puig'))
    print(sim_taxicab(critics, 'Toby', 'Claudia Puig'))

    print "Lisa Rose Top Matches"
    print(topMatches(critics, 'Lisa Rose', n = 3))
    print(topMatches(critics, 'Lisa Rose', n = 3, similarity = sim_distance))
    print(topMatches(critics, 'Lisa Rose', n = 3, similarity = sim_cosine))
    print(topMatches(critics, 'Lisa Rose', n = 3, similarity = sim_jaccard))
    print(topMatches(critics, 'Lisa Rose', n = 3, similarity = sim_tanimoto))
    print(topMatches(critics, 'Lisa Rose', n = 3, similarity = sim_taxicab))
    print "Toby Top Matches"
    print(topMatches(critics, 'Toby', n = 3))
    print(topMatches(critics, 'Toby', n = 3, similarity = sim_distance))
    print(topMatches(critics, 'Toby', n = 3, similarity = sim_cosine))
    print(topMatches(critics, 'Toby', n = 3, similarity = sim_jaccard))
    print(topMatches(critics, 'Toby', n = 3, similarity = sim_tanimoto))
    print(topMatches(critics, 'Toby', n = 3, similarity = sim_taxicab))

    print '----Recommendations for Michael Phillips'
    print(getRecommendations(critics, 'Michael Phillips'))
    print(getRecommendations(critics, 'Michael Phillips', similarity = sim_distance))
    print(getRecommendations(critics, 'Michael Phillips', similarity = sim_cosine))
    print(getRecommendations(critics, 'Michael Phillips', similarity = sim_jaccard))
    print(getRecommendations(critics, 'Michael Phillips', similarity = sim_tanimoto))
    print(getRecommendations(critics, 'Michael Phillips', similarity = sim_taxicab))
    print '----Recommendations for Toby'
    print(getRecommendations(critics, 'Toby'))
    print(getRecommendations(critics, 'Toby', similarity = sim_distance))
    print(getRecommendations(critics, 'Toby', similarity = sim_cosine))
    print(getRecommendations(critics, 'Toby', similarity = sim_jaccard))
    print(getRecommendations(critics, 'Toby', similarity = sim_tanimoto))
    print(getRecommendations(critics, 'Toby', similarity = sim_taxicab))

    movies = transformPrefs(critics)

    print 'Superman Returns Top Matches'
    print(topMatches(movies, 'Superman Returns'))
    print(topMatches(movies, 'Superman Returns', similarity = sim_distance))
    print(topMatches(movies, 'Superman Returns', similarity = sim_cosine))
    print(topMatches(movies, 'Superman Returns', similarity = sim_jaccard))
    print(topMatches(movies, 'Superman Returns', similarity = sim_tanimoto))
    print(topMatches(movies, 'Superman Returns', similarity = sim_taxicab))
    print 'Just My Luck Top Matches'
    print(topMatches(movies, 'Just My Luck'))
    print(topMatches(movies, 'Just My Luck', similarity = sim_distance))
    print(topMatches(movies, 'Just My Luck', similarity = sim_cosine))
    print(topMatches(movies, 'Just My Luck', similarity = sim_jaccard))
    print(topMatches(movies, 'Just My Luck', similarity = sim_tanimoto))
    print(topMatches(movies, 'Just My Luck', similarity = sim_taxicab))

    '''
    print(getRecommendations(movies, 'Superman Returns'))
    print(getRecommendations(movies, 'Just My Luck'))
    print(getRecommendations(movies, 'Superman Returns', similarity = sim_distance))
    print(getRecommendations(movies, 'Just My Luck', similarity = sim_distance))
    print(getRecommendations(movies, 'Superman Returns', similarity = sim_taxicab))
    print(getRecommendations(movies, 'Just My Luck', similarity = sim_taxicab))
    '''

    print '====='
    print '---Similarity Pearson---'
    itemsim = calculateSimilarItems(critics)
    #print itemsim
    print getRecommendedItems(critics, itemsim, 'Michael Phillips')
    print getRecommendedItems(critics, itemsim, 'Toby')
    print '---Similarity Distance---'
    itemsim = calculateSimilarItems(critics, similarity = sim_distance)
    #print itemsim
    print getRecommendedItems(critics, itemsim, 'Michael Phillips')
    print getRecommendedItems(critics, itemsim, 'Toby')
    print '---Similarity Cosine---'
    itemsim = calculateSimilarItems(critics, similarity = sim_cosine)
    #print itemsim
    print getRecommendedItems(critics, itemsim, 'Michael Phillips')
    print getRecommendedItems(critics, itemsim, 'Toby')
    print '---Similarity jaccard---'
    itemsim = calculateSimilarItems(critics, similarity = sim_jaccard)
    #print itemsim
    print getRecommendedItems(critics, itemsim, 'Michael Phillips')
    print getRecommendedItems(critics, itemsim, 'Toby')
    print '---Similarity tanimoto---'
    itemsim = calculateSimilarItems(critics, similarity = sim_tanimoto)
    #print itemsim
    print getRecommendedItems(critics, itemsim, 'Michael Phillips')
    print getRecommendedItems(critics, itemsim, 'Toby')
    print '---Similarity taxicab---'
    itemsim = calculateSimilarItems(critics, similarity = sim_taxicab)
    #print itemsim
    print getRecommendedItems(critics, itemsim, 'Michael Phillips')
    print getRecommendedItems(critics, itemsim, 'Toby')

