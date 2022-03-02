# Mutation-analysis
This is intended to fill out data that was mined from the dmd.nl database on Muscular Dystrophy with relevant physical characteristics and assess any trends in physical characteristics and the likelihood of being a disease causing mutation.


The data was mined from the dmd.nl database with necessary permissions. Also during the mining of this data, I comitted a couple of updates to the repository to improve the accuracy of the documented mutations (through comparing the amino acid 1 letter, 3 letter,  DNA sequence predictions, and PDB sequences).

Then data concerning physical characteristics such as amino acid size, polarity, acidity/basicity are all filled into the table based on empirical models that are cited within the code including amino acide heat capacity.

Exploratory analysis through pariwise graphing and correlation are all followed by significance testing through ANOVA and t-tests
comparing disease causing and non-disease causing mutations to inform studies on the Dystrophin protein for Biophysical purposes.


# Ultimately, this analysis conveyed a significant difference in the predicted change in heat capacity
between disease causing and non-disease causing mutations that informed our lab to use Differential Scanning calorimetry to study this protein as it is the only method that can detect a change in heat capacity directly.
