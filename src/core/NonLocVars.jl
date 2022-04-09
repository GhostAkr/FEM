"""
    NonLocPars

Parameters of non-local simulation model. Contains following fields:
1. `locimpact_coeff::Number`: coefficient denoting local impact;
2. `nonlocimpact_coeff::Number`: coefficient denoting non-local impact;
3. `impactdist::Number`: impact distance of non-local effects.
"""
struct NonLocPars
    locimpact_coeff::Number
    nonlocimpact_coeff::Number
    impactdist::Number
end