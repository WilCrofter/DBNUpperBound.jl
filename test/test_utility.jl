
#= using code from KM's Python 2.7 repository:
https://github.com/km-git-acc/dbn_upper_bound.
The import statement executes the following
#check phi_decay values
print phi_decay(0.001)
print phi_decay(0.01)
print phi_decay(0.1)
print phi_decay(0.5)
print phi_decay(1)

Python 2.7.14 (default, Sep 20 2017, 01:25:59) 
[GCC 7.2.0] on linux2
Type "help", "copyright", "credits" or "license" for more information.
>>> import utility
0.446680170237
0.445026555345
0.304274852667
1.37781394064e-07
5.10200133902e-70
>>> 
=#

@test phi_decay(0.001) ≈ 0.446680170237
@test phi_decay(0.01) ≈ 0.445026555345
@test phi_decay(0.1) ≈ 0.304274852667
@test phi_decay(0.5) ≈ 1.37781394064e-07
@test phi_decay(1) ≈ 5.10200133902e-70



