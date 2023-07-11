from django.shortcuts import render
from django.http import HttpResponse

import math
import scipy.stats as ss
 
def d(S0, sigma, r, T, K):
    return (math.log(S0 / K) + (r + sigma ** 2 / 2) * T) / (sigma * math.sqrt(T))

def  euro(S0, sigma, r, T, K, option_type):
    if option_type == "C":
        return round(S0 * ss.norm.cdf(d(S0, sigma, r, T, K)) - K * math.exp(-r * T) * ss.norm.cdf(d(S0, sigma, r, T, K)-sigma*math.sqrt(T)), 3)
    else:
        return round((-S0) * ss.norm.cdf((-1)*d(S0, sigma, r, T, K)) + K * math.exp(-r * T) * ss.norm.cdf((-1)*d(S0, sigma, r, T, K)+sigma*math.sqrt(T)), 3)


def index(request):
    return render(request, 'Blackscholes.html')


def blackScholes(request, S, sigma, r, T, K, call_put):
    c =  euro(int(S), float(sigma), float(r), float(T), float(K), str(call_put))
    r = HttpResponse(str(c))
    return r
