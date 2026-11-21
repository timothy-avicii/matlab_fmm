# -*- coding: utf-8 -*-
"""
Created on Sun Jan 26 10:09:57 2025
modules also require pandas matplotlib 
@author: atim9
"""
import numpy as np
import torch
import pybbfmm
import aljpy

# r = 來源點r′和目標點r之間的距離
# alpha 前向散射距離
# beta  後向散射距離
# eta   後向散射與前向散射能量比值
# 雙高斯函數
def double_gaussian(sources, targets):
    r = torch.norm(sources - targets, dim=-1) 
    alpha = 10
    beta = 500
    eta = 0.7
    return (1 / (np.pi * (1 + eta))) * (
        (1 / alpha**2) * torch.exp(-r**2 / alpha**2) +
        (eta / beta**2) * torch.exp(-r**2 / beta**2)
    )

# gamma 中間散射區域距離
# eta_prime  中間散射與前向散射能量比值
# 三高斯函數
def triple_gaussian(sources, targets):
    r = torch.norm(sources - targets, dim=-1)   
    alpha = 50
    beta = 1000
    gamma = 500
    eta = 0.3
    eta_prime = 0.15
    return (1 / (np.pi * (1 + eta + eta_prime))) * (
        (1 / alpha**2) * np.exp(-r**2 / alpha**2) +
        (eta / beta**2) * np.exp(-r**2 / beta**2) +
        (eta_prime / gamma**2) * np.exp(-r**2 / gamma**2)
    )

# 雙高斯函數 + 指數項
# eta_double_prime遠場尾部效應
# gamma_prime 遠場尾部效應距離
def double_gaussian_plus_exponential(sources, targets):
    r = torch.norm(sources - targets, dim=-1) 
    alpha = 50
    beta = 1000
    gamma = 500
    eta = 0.3
    eta_prime = 0.1
    return (1 / (np.pi * (1 + eta + eta_prime))) * (
        (1 / alpha**2) * torch.exp(-r**2 / alpha**2) +
        (eta / beta**2) * torch.exp(-r**2 / beta**2) +
        (eta_prime / (2 * gamma**2)) * torch.exp(-r / gamma)
    )

# 三高斯函數 + 指數項
def triple_gaussian_plus_exponential(sources, targets):
    r = torch.norm(sources - targets, dim=-1)
    alpha = 50
    beta = 1000
    gamma = 500  
    gamma_prime = 1500  
    eta = 0.3
    eta_prime = 0.15
    eta_double_prime = 0.05
    return (1 / (np.pi * (1 + eta + eta_prime + eta_double_prime))) * (
        (1 / alpha**2) * torch.exp(-r**2 / alpha**2) +
        (eta / beta**2) * torch.exp(-r**2 / beta**2) +
        (eta_prime / gamma**2) * torch.exp(-r**2 / gamma**2) +
        (eta_double_prime / (2 * gamma_prime**2)) * torch.exp(-r / gamma_prime)
    )

#chebyshev degree不要太大
def matlab_fmm(layout, cheb_degree = 4, Gaussian = 0):
    source_points = torch.nonzero(layout, as_tuple=False) #找到所有曝光點的位置 
    dose = torch.ones(source_points.size(0))
     
    #葉節點最大曝光點數在其他程式內計算
    rows = torch.arange(0, layout.size(0))
    cols = torch.arange(0, layout.size(1))
    target_points = torch.cartesian_prod(rows, cols)
    
    input_data = aljpy.dotdict({
        'sources': source_points,
        'targets': target_points, 
        'charges' : dose, 
    })
    if Gaussian == 0 :
        input_data['kernel'] = lambda a, b: double_gaussian(a, b)
    elif Gaussian == 1 :
        input_data['kernel'] = lambda a, b: triple_gaussian(a, b)
    elif Gaussian == 2 :
        input_data['kernel'] = lambda a, b: double_gaussian_plus_exponential(a, b)
    elif Gaussian == 3 :
        input_data['kernel'] = lambda a, b: triple_gaussian_plus_exponential(a, b)
    else:
        print('invalid Gaussian kernel')
    #Gaussian_kernel    
    eva = pybbfmm.solve(input_data, cheb_degree)
    
    result_matrix = eva.view(layout.size(0),layout.size(1))
    return result_matrix




