#!/usr/bin/env python
# coding: utf-8

def produtos_naturais (lotus, anotacao):
    import pandas as pd

    df = lotus.merge(anotacao, on=['inchi'], how='outer', suffixes=['', '_'], indicator=True)
    both = df.loc[df['_merge']== 'both']
    return both