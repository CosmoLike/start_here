import numpy as np

#Fixed seed so we always get the same blinding
seed=1062016

#Numpy RNG; should be repeatable
rng=np.random.RandomState(seed)

#Ensure that the RNG has not changed with
#different versions of numpy or different
#machines. If this test value is different
#than the value below, which I just ran on 
#my machine, then we have a problem
#from the one on my machine then we have a problem
test_value = rng.uniform()
assert abs(0.0541095275873-test_value)<1e-9

#The parameters we want to blind and the magnitude of 
#the shifts we will apply to them
blinding_scales = [
    ("omegam", 0.1),
    ("sigma8", 0.2),
    ("h0", 0.2),
    ("ns", 0.1),
    ("w0", 0.4),
    ("omegab", 0.05),
]

def normalize_name(name):
    return name.replace("_", "").lower()

blinded_parameters = [b[0] for b in blinding_scales]

#generate blinding factors, different for each parameter,
#between -1 and 1 uniformly.
#Please do not look at these factors
blinding_factors = rng.uniform(size=len(blinding_scales))*2-1


def blind_parameters(names, values):
    output = []
    assert len(names)==len(values)
    for name, value in zip(names, values):
        try:
            i = blinded_parameters.index(normalize_name(name))
        except ValueError:
            i = -1
        if i<0:
            output.append(value)
        else:
            f = blinding_factors[i]
            s = blinding_scales[i][1]
            output.append(value + f*s)
    return np.array(output)




def unblind_parameters(names, values):
    output = []
    assert len(names)==len(values)
    for name, value in zip(names, values):
        try:
            i = blinded_parameters.index(normalize_name(name))
        except ValueError:
            i = -1
        if i<0:
            output.append(value)
        else:
            f = blinding_factors[i]
            s = blinding_scales[i][1]
            output.append(value - f*s)
    return np.array(output)



if __name__ == '__main__':
    names = ['omega_m', 'sigma_8', 'shear_m_0', 'bias_1', 'w0']
    values = [0.3, 0.8, 0.0, 1.5, -1.0]
    values2 = blind_parameters(names, values)
    values3 = unblind_parameters(names, values2)

    print values
    print values2
    print values3