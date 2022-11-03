
The list of operations was extracted from https://numpy.org/doc/stable/user/basics.indexing.html


The following are not (yet) working:

```
x[:, np.newaxis, :, :].shape
x[:, None, :, :].shape
x[:, np.newaxis] + x[np.newaxis, :]
x[~np.isnan(x)]
x[x < 0]
x[x < 0] += 20
x[2:7] = 1
x[2:7] = np.arange(5)
x[np.array([1, 1, 3, 1])] += 1
```