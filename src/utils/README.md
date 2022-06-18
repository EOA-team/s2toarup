# rtm_inv

`rtm_inv` is a light-weight module to perform lookup-table based inversion
of the radiative transfer models ProSAIL and (in the future) SPART.

`rtm_inv` is undergoing current developments, therefore, it might be not stable.

The sub-package provided here is for reproducible reasons but might differ from
any future versions of `rtm_inv`.

`rtm_inv` requires the following packages:

```
pandas
numpy
matplotlib
prosail
spectral
numba
```

You might have to explicitely install `prosail` and `spectral` if you want to create LUTs on your own.
