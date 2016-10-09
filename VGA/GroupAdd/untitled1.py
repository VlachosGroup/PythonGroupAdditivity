# -*- coding: utf-8 -*-
"""
Created on Fri Oct 07 14:04:25 2016

@author: Gu
"""
import os
import yaml
path = 'benson'
__file__
base_path = os.path.join(os.path.split(os.path.split(os.path.split(__file__)[0])[0])[0],'data')
path = os.path.join(base_path,path)
path = os.path.join(path,'library.yaml')
with open(path) as f:
    library = yaml.load(f)