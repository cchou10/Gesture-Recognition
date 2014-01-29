ECE5470
=======

Final project from ECE5470: Computer Vision. Vision system to read a small subset of American Sign Language.

Abstract:
This project approaches the classication aspect of gesture recognition. It recognizes numbers 1
through 10 in American Sign Language (ASL). The general algorithm rst thresholds the greyscale
image. It then calculates center of intensity (COI) in order to predict locations for ngertips. At
these locations it draws bounding boxes as a mask and counts how many hand pixels are in each.
It can then determine the status of the ngers and thus recognize the ASL number. However, due
to problems with lighting and the general dataset other considerations such as thresholding and
limiting our region of interest wet taken. This gesture recognition has many applications such as an
interactive command communication with robots and electronic devices, security passwords, and
the obvious application for translating ASL.

