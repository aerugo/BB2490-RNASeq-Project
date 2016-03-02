I've moved my diary to the project wiki, https://rawgit.com/aerugo/BB2490-RNASeq-Project/wiki/Tobias'-Logbook

# Diary - Tobias Frick

## Week one

### Captains log, stardate 1502.17 - Startup

Today we planned some meetings for the coming weeks and later, Hugi introduced me and Amy to PyCharm ans GitHub.
We set up a working directory and also tried to push/pull to the working directory on GitHub. I had some minor problems
with setting the rig up since I'm using Windows but with a little bit of work it got sorted out. We (aka Hugi) also 
checked that GitHub is compatible with UppMax, which it seemed to be. As a last thing to mention, I briefed the group I will be busy taking a lab course the next week which will be occupying most of my time unfortunately.

### Captains log, stardate 1502.18 - Project Description

Today I've created a facebook group for communicating with the group. More importantly the project discription and (preliminary) title arrived today so I created a draft of the presentation using the project description and google drive. I also made it clear to the group that they are free to edit it if they wish. PS: Ideally the group would have created the presentation from scratch together, however lack of time before the presentation in combination with poorly matching schedules did not allow this according to me.

### Captains log, stardate 1502.19 - Project Presentation

Today we've presented our projects to the other groups and also listened to the other group presentations about their project. After that rather short session we had a meeting with Olof (our supervisor) and discussed the project description a little more. As a last thing for the group today, we talked a little about what I can do this weekend to make up fo my upcoming lab week. It was decided I should look into studies regarding what softwares to use ans also I should try to find out more about the hypothosis about reaching sequencing saturation.

Supplemental: I did some minor searching on PubMed to try and find out more about the hypothesis of depths saturation but did not find anything relevant. 

### Captains log, stardate 1502.21 - Literature study

I have been searching and reading articles regarding the project, mainly focused on finding comparative studies for software (both for DE analysis and allelic bias). Unfortunately I still have not found any documentation of the saturation hypothesis. This is what reasonable studies I've found so far:

**DE software comparison studies:**
doi: 10.1186/1471-2105-14-91
doi: 10.1093/bib/bbt086
doi: 10.1371/journal.pone.0103207
doi:10.3390/bios3030238
doi: 10.3732/ajb.1100340

**ASE software articles:**
doi: 10.1038/nmeth.3582
doi: 10.1093/bioinformatics/btu802
doi: 10.1186/s13059-014-0405-3
doi: 10.1038/srep21134

## Week two

**NB:** I've been busy cultivating cells this week so I havn't had time to work on the project and therefore worked the weeekend before and will also put in some extra work the coming weekend.

### Captains log, stardate 1502.26 - Seminar, project discussion and initial runs

The first thing today was attending to the seminar and also presenting our progress since last friday. During this time, Sailendra also briefed me on what he and Hugi did last tuesday to make sure I did not miss anything important when reading their results/diaries. We had some minor discussion regarding our problems with Johannes Alneberg.

After lunch Sailendra and I met up to work on the project. We started by discussing what questions we actually wanted to answer with this project, mainly focused around the random sampling and the problems of the unevenly distributed reads in different flow cells. Somewhere aruond here, Hugi arrived and he joined the discussion. We ended at the point where we were unsure of whether we were supposed to run the analysis on the **pooled genome subsamplings**, **individual subsampling tests** or if we should try to make all the reads have an **even distribution** by **discarding data**. We decided that we should talk to Olof regarding our problems regarding which question to answer and the "strangeness" of the unevenly distributed reads in the flowcells. We also asked about the hypothesis about saturation in hopes of being provided with some articles on the subject.

https://rawgit.com/aerugo/BB2490-RNASeq-Project/master/results/Readcounts.png
**[INPUT GRAPH OF READ DISTRIBUTION]**

When the email to Olof was sent we started working. Sailendra and Hugi started on the mapping and also looked a little bit at running STAR whilst I looked further into the DE software analysis. After some time of reading I recommended DESeq as a suitable tool for our analysis since it had been recommended by studies and seemed to perform very well. Granted, there was several suitable software but there were only a few which performed at the top level and of these few top-performers (**DESeq**, **edgeR** and **limma**) DESeq was the most commonly mentioned. Being the most common software compared might not be a merit in itself but it also showed to be more documented and therefore I thought that to be a good tiebraker as for now. Furthermore, edgeR was shown to be more liberal than DESeq but also accepting more FP and having an increased FPR (FP rate) than DESeq which is according to my opinion something that shows it is not as good as DESeq (but it should be noted that is just my personal preference but the group seemed to agree).

**[LINK TO THE TWO MAIN STUDIES IN RESOURCES]**

One of the articles also mentioned some things about sequencing depths and more specifically what they called **unbalanced sequencing depths**, which seemed to be what we have in our situtaion. In this study they showed that DESeq (and other DE-software) handled in group variations which would mean that potentially our in group variation will not affect the DE results, which would be nice. It should be noted that a **between group variation** seemed to impact the DE analysis negatively which would mean our current 2:nd sample  (high sequencing depths) sample would be impaired in its ability to find DE genes (DEG).

The last thing I did was to star looking at how GeneiASE works, what it does and what makes it interesting. So far I've gotten some information like: it **compensates for mapping bias** by using a modified binominal distribution, it analyses the icd-ASE which means it cannot only see which allele is exressed more but also make a comparison between two individual conditions. So basically it produces the allelic inbalance for a specific alles whilst also seing how that balance changes between conditions. By now I have to start to dig into the supplemental information in order to find out more about the method GeneIASE uses to counteract the mapping bias and hopefully we can have a discussion about it at our next meeting with Olof.
