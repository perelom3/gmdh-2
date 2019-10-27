# Read data
# Split to train test exam parts

# --- PREVIEW ----
# Data
# FFS over data

# ---- INPUT FROM USER -----
# frequencies amount(m) and interval
# stop threshold 
# GA operators
# start amount of applicants (M)
# k - interations to double check
# F - freedom param
# mutation probability
# cross over probability


# ---- GMDH --------

# Generate M applicants with m chromosomes
# for each applicant
# -- compute least squares method(LSM)
# -- compute external criterion(EC) for model on test data
# select H < M best applicants for next breed by EC
# Save F models
# breed new M > H applacants from H selected
# Check if GA should stop by absolute value of difference between current best EC and previous, but keep going for k steps to double check

# ---- Result ------
# print best 5 models
# plot best model and data
# calc criterions on exam data