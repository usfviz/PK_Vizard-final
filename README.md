 # Team
 Phil Rotella
 
 Brad Kenstler
 
 # Packages
* `dplyr`
* `googleVis`
* `ggvis`

 # Dataset
 Our dataset for this application comes from the U.S. Food & Drug Administrations open API.
 The data in question is a report of "Analyzed Events" (AE), which are documented complications patients have experienced in relation to some drug.
 The raw events have been binned into 1000 AE categories, such as "Dizziness" or "Loss of libido."
 
 We've taken this data and organized it to various frequencies, including frequency of events for a specific AE or frequency of any event for a specific drug.
 From this we've been able to put together the following application, which looks at the relationships between these two quantities.
 
 # Interface
 
 Our application is a 4-panel visualization of this data.
 
<img src="images/shiny1.png"/>

The first panel shows the top ten AE's by frequency for the specified drug.

<img src="images/shiny2.png"/>

The second panel shows the top ten Drugs by frequency for a given AE.

<img src="images/shiny3.png"/>

The third panel is a radial network. It contains the ten drugs with the most recorded events, and for each of those drugs, the ten most frequent AE's

<img src="images/shiny6.png"/>

The third panel shows a table of AE's for each drug, ranked by the LLR. The LLR is the log-likelihood ratio statistic for that given drug and AE combination. The LLR is a statistic which expresses how many times more likely the observed data is under one model than the other. In our context, we're using it as a measure of how different the AE distribution for a given drug is compared to all drugs; in other words, do these events indicate that this drug is a plausible cause for high occurrences of an AE or is it just a "normal" occurrence.

The quantity we'd like to have to confirm this is a p-value. The fina version of this project may include bootstrapped p-values for each combination.

<img src="images/shiny5.png"/>

Our last panel is a heatmap of drug/AE combinations for a subset of drugs. The interface allows you to select which drug to sort the AE's on the y-axis by.

