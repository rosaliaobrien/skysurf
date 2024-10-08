# SKYSURF Completeness GUI

To use the SKYSURF Completeness GUI, simply run

```
python completeness_calculator_gui.py
```

Specify the HST camera, filter, object size (in arcsec), background level (in counts per sec), and exposure time (in sec).
A good example is using size of 1 arcsec, 1000 seconds of exposure time, and 0.1 for background. The 50% completeness magnitude is then displayed.

Requirements: PyQt5 and sys packages
