from django.urls import include, re_path
from django.contrib import admin
admin.autodiscover()

from django.shortcuts import render
from django.http import HttpResponse

from functions.Blackscholes_views import blackScholes, index

import math
import scipy.stats as ss
from math import *


 
def binomial_tree(S0, K, T, r, sigma, N, call_put):
    t = T/N
    down = math.exp(-sigma * math.sqrt(t))
    up = math.exp(sigma * math.sqrt(t))
    p = (math.exp(r*t) - down)/(up - down)
    c = 0
    for i in range(N+1):
        
        node_prob = comb(N, i) * (p)**i * (1-p)**(N-i)
        
        ST = S0 * (up)**i * (down)**(N-i)
        if call_put == 'P':
            c += max(K-ST, 0) * node_prob
        elif call_put == 'C':
            c += max(ST-K, 0) * node_prob
        else:
            raise ValueError("Must be 'C' or 'P'")
    return round(c * math.exp(-r*T),3)


def binomialIndex(request):
    return render(request, 'BinomialModel.html')


def Binomial(request, S, K, T, r, sigma, N, call_put):
    c = binomial_tree(float(S), float(K), float(T), float(r), float(sigma), int(N), str(call_put))
    r = HttpResponse(str(c))
    return r

urlpatterns = [
    #Black-Scholes
    re_path(r'^$', index, name='index'),
    re_path(r'^catalog', index, name='index'),
    re_path(r'^blackScholes_index', index, name='index'),
    re_path(r'^blackScholes/S(.+)sigma(.+)r(.+)T(.+)K(.+)call_put(.+)/$', blackScholes, name='blackScholes'),

    #Binomial model
    re_path(r'^binomial_index', binomialIndex, name='binomial'),
    re_path(r'^binomial/S(.+)K(.+)T(.+)r(.+)sigma(.+)N(.+)call_put(.+)/$', Binomial, name='n-step Binomial'),
]
