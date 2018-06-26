pro pvalue_binomial

; Ts >= h
h_3 = 8.88d
p_h_3 = 2.69e-03

N_slide = 6.
k_slide = 2.

Prob = binomial(k_slide, N_slide, p_h_3)

print, 'Using the BINOMIAL tool:'
print, 'Probability =', Prob

quantile = (inverf(1. - Prob))*sqrt(2.)

print, 'Sigma = ', quantile

end
