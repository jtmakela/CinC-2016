# Time and Frequency -Based Approach to Heart Sound Segmentation and Classification

This is our (Jarno Mäkelä, Heikki Väänänen) official entry for Classification of Normal/Abnormal Heart Sound Recordings: [the PhysioNet / Computing in Cardiology Challenge 2016](https://www.physionet.org/challenge/2016/).

There is also our final full paper as presented in [Computing in Cardiology 2016, Volume 43](http://www.cinc.org/archives/2016/index.html#session6-6)

In this study, we propose a decision tree classifier of heart sound signals. We determined repetitive fundamental heart sound segments based on adaptive similarity value clusterization of the sound signal, and we created a set of filters for decision tree parametrization. Using the filters together with inter-segment timings, we created three sets of markers: a set utilizing both S1 and S2 identification, a set where only one segment was identified, and a set without any identified segment. An individual classification tree was trained for each marker set.

As a result, our classifier attained sensitivity (Se) of 0.66 and specificity (Sp) of 0.92 and overall score of 0.79 for a hidden random (revised) subset.

![alt tag](https://github.com/jtmakela/CinC-2016/blob/master/s1.average.png)
