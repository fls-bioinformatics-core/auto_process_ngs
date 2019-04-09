#     mockqcdata.py: module providing mock Illumina data for testing
#     Copyright (C) University of Manchester 2012-2019 Peter Briggs
#
########################################################################

FASTQC_0_11_3 = {
    'fastqc_data.txt': """##FastQC	0.11.3
>>Basic Statistics	pass
#Measure	Value
Filename	%(fastq)s
File type	Conventional base calls
Encoding	Sanger / Illumina 1.9
Total Sequences	10
Sequences flagged as poor quality	0
Sequence length	250
%%GC	51
>>END_MODULE
>>Per base sequence quality	fail
#Base	Mean	Median	Lower Quartile	Upper Quartile	10th Percentile	90th Percentile
1	27.6	0.0	0.0	0.0	0.0	0.0
2	25.6	0.0	0.0	0.0	0.0	0.0
3	22.0	0.0	0.0	0.0	0.0	0.0
4	26.8	0.0	0.0	0.0	0.0	0.0
5	28.4	0.0	0.0	0.0	0.0	0.0
6	31.3	0.0	0.0	0.0	0.0	0.0
7	29.4	0.0	0.0	0.0	0.0	0.0
8	29.7	0.0	0.0	0.0	0.0	0.0
9	31.9	0.0	0.0	0.0	0.0	0.0
10-14	31.639999999999997	0.0	0.0	0.0	0.0	0.0
15-19	34.48	0.0	0.0	0.0	0.0	0.0
20-24	30.339999999999996	0.0	0.0	0.0	0.0	0.0
25-29	34.019999999999996	0.0	0.0	0.0	0.0	0.0
30-34	34.839999999999996	0.0	0.0	0.0	0.0	0.0
35-39	35.6	0.0	0.0	0.0	0.0	0.0
40-44	34.56	0.0	0.0	0.0	0.0	0.0
45-49	36.36	0.0	0.0	0.0	0.0	0.0
50-54	35.279999999999994	0.0	0.0	0.0	0.0	0.0
55-59	35.580000000000005	0.0	0.0	0.0	0.0	0.0
60-64	35.56	0.0	0.0	0.0	0.0	0.0
65-69	34.58	0.0	0.0	0.0	0.0	0.0
70-74	34.67999999999999	0.0	0.0	0.0	0.0	0.0
75-79	34.17999999999999	0.0	0.0	0.0	0.0	0.0
80-84	32.220000000000006	0.0	0.0	0.0	0.0	0.0
85-89	34.22	0.0	0.0	0.0	0.0	0.0
90-94	34.49999999999999	0.0	0.0	0.0	0.0	0.0
95-99	34.3	0.0	0.0	0.0	0.0	0.0
100-104	33.52	0.0	0.0	0.0	0.0	0.0
105-109	32.06	0.0	0.0	0.0	0.0	0.0
110-114	28.2	0.0	0.0	0.0	0.0	0.0
115-119	30.639999999999997	0.0	0.0	0.0	0.0	0.0
120-124	28.840000000000003	0.0	0.0	0.0	0.0	0.0
125-129	28.96	0.0	0.0	0.0	0.0	0.0
130-134	25.860000000000003	0.0	0.0	0.0	0.0	0.0
135-139	27.660000000000004	0.0	0.0	0.0	0.0	0.0
140-144	26.78	0.0	0.0	0.0	0.0	0.0
145-149	23.96	0.0	0.0	0.0	0.0	0.0
150-154	23.9	0.0	0.0	0.0	0.0	0.0
155-159	25.2	0.0	0.0	0.0	0.0	0.0
160-164	22.020000000000003	0.0	0.0	0.0	0.0	0.0
165-169	21.639999999999997	0.0	0.0	0.0	0.0	0.0
170-174	17.72	0.0	0.0	0.0	0.0	0.0
175-179	16.36	0.0	0.0	0.0	0.0	0.0
180-184	15.540000000000001	0.0	0.0	0.0	0.0	0.0
185-189	15.88	0.0	0.0	0.0	0.0	0.0
190-194	19.82	0.0	0.0	0.0	0.0	0.0
195-199	16.419999999999998	0.0	0.0	0.0	0.0	0.0
200-204	17.76	0.0	0.0	0.0	0.0	0.0
205-209	15.319999999999999	0.0	0.0	0.0	0.0	0.0
210-214	16.9	0.0	0.0	0.0	0.0	0.0
215-219	17.080000000000002	0.0	0.0	0.0	0.0	0.0
220-224	14.839999999999998	0.0	0.0	0.0	0.0	0.0
225-229	11.760000000000002	0.0	0.0	0.0	0.0	0.0
230-234	13.74	0.0	0.0	0.0	0.0	0.0
235-239	12.54	0.0	0.0	0.0	0.0	0.0
240-244	8.14	0.0	0.0	0.0	0.0	0.0
245-249	3.7399999999999998	0.0	0.0	0.0	0.0	0.0
250	4.5	0.0	0.0	0.0	0.0	0.0
>>END_MODULE
>>Per tile sequence quality	pass
#Tile	Base	Mean
1101	1	0.0
1101	2	0.0
1101	3	0.0
1101	4	0.0
1101	5	0.0
1101	6	0.0
1101	7	0.0
1101	8	0.0
1101	9	0.0
1101	10-14	0.0
1101	15-19	0.0
1101	20-24	0.0
1101	25-29	0.0
1101	30-34	0.0
1101	35-39	0.0
1101	40-44	0.0
1101	45-49	0.0
1101	50-54	0.0
1101	55-59	0.0
1101	60-64	0.0
1101	65-69	0.0
1101	70-74	0.0
1101	75-79	0.0
1101	80-84	0.0
1101	85-89	0.0
1101	90-94	0.0
1101	95-99	0.0
1101	100-104	0.0
1101	105-109	0.0
1101	110-114	0.0
1101	115-119	0.0
1101	120-124	0.0
1101	125-129	0.0
1101	130-134	0.0
1101	135-139	0.0
1101	140-144	0.0
1101	145-149	0.0
1101	150-154	0.0
1101	155-159	0.0
1101	160-164	0.0
1101	165-169	0.0
1101	170-174	0.0
1101	175-179	0.0
1101	180-184	0.0
1101	185-189	0.0
1101	190-194	0.0
1101	195-199	0.0
1101	200-204	0.0
1101	205-209	0.0
1101	210-214	0.0
1101	215-219	0.0
1101	220-224	0.0
1101	225-229	0.0
1101	230-234	0.0
1101	235-239	0.0
1101	240-244	0.0
1101	245-249	0.0
1101	250	0.0
>>END_MODULE
>>Per sequence quality scores	pass
#Quality	Count
12	1.0
13	0.0
14	0.0
15	0.0
16	1.0
17	0.0
18	1.0
19	0.0
20	0.0
21	0.0
22	1.0
23	0.0
24	0.0
25	0.0
26	0.0
27	1.0
28	0.0
29	1.0
30	1.0
31	1.0
32	0.0
33	2.0
>>END_MODULE
>>Per base sequence content	fail
#Base	G	A	T	C
1	40.0	30.0	10.0	20.0
2	20.0	30.0	40.0	10.0
3	20.0	40.0	40.0	0.0
4	10.0	20.0	20.0	50.0
5	0.0	60.0	30.0	10.0
6	40.0	20.0	10.0	30.0
7	10.0	20.0	40.0	30.0
8	10.0	20.0	30.0	40.0
9	40.0	10.0	0.0	50.0
10-14	26.0	22.0	20.0	32.0
15-19	24.0	46.0	12.0	18.0
20-24	21.73913043478261	26.08695652173913	21.73913043478261	30.434782608695656
25-29	20.0	34.0	18.0	28.000000000000004
30-34	22.0	14.000000000000002	28.000000000000004	36.0
35-39	30.0	34.0	18.0	18.0
40-44	22.0	22.0	30.0	26.0
45-49	28.000000000000004	24.0	24.0	24.0
50-54	34.0	20.0	20.0	26.0
55-59	24.0	24.0	12.0	40.0
60-64	12.0	28.000000000000004	32.0	28.000000000000004
65-69	28.000000000000004	16.0	34.0	22.0
70-74	22.0	26.0	24.0	28.000000000000004
75-79	32.0	16.0	20.0	32.0
80-84	24.0	22.0	22.0	32.0
85-89	14.000000000000002	24.0	34.0	28.000000000000004
90-94	28.000000000000004	28.000000000000004	28.000000000000004	16.0
95-99	28.000000000000004	26.0	10.0	36.0
100-104	20.0	16.0	34.0	30.0
105-109	18.0	18.0	24.0	40.0
110-114	30.0	20.0	18.0	32.0
115-119	18.0	28.000000000000004	30.0	24.0
120-124	26.0	28.000000000000004	18.0	28.000000000000004
125-129	14.000000000000002	26.0	22.0	38.0
130-134	32.0	30.0	22.0	16.0
135-139	22.0	22.0	32.0	24.0
140-144	20.0	36.0	16.0	28.000000000000004
145-149	32.0	18.0	28.000000000000004	22.0
150-154	26.0	26.0	28.000000000000004	20.0
155-159	16.0	22.0	24.0	38.0
160-164	26.0	20.0	36.0	18.0
165-169	30.0	24.0	18.0	28.000000000000004
170-174	28.000000000000004	16.0	24.0	32.0
175-179	24.0	34.0	14.000000000000002	28.000000000000004
180-184	24.0	24.0	20.0	32.0
185-189	22.0	24.0	24.0	30.0
190-194	16.0	14.000000000000002	36.0	34.0
195-199	22.0	30.0	20.0	28.000000000000004
200-204	24.0	30.0	26.0	20.0
205-209	24.0	22.0	26.0	28.000000000000004
210-214	20.0	22.0	26.0	32.0
215-219	18.0	24.0	28.000000000000004	30.0
220-224	18.0	24.0	22.0	36.0
225-229	26.0	24.0	20.0	30.0
230-234	20.0	30.0	18.0	32.0
235-239	28.000000000000004	16.0	18.0	38.0
240-244	28.000000000000004	32.0	22.0	18.0
245-249	22.0	34.0	12.0	32.0
250	30.0	40.0	10.0	20.0
>>END_MODULE
>>Per sequence GC content	fail
#GC Content	Count
0	0.0
1	0.0
2	0.0
3	0.0
4	0.0
5	0.0
6	0.0
7	0.0
8	0.0
9	0.0
10	0.0
11	0.0
12	0.0
13	0.0
14	0.0
15	0.0
16	0.0
17	0.0
18	0.0
19	0.0
20	0.0
21	0.0
22	0.0
23	0.0
24	0.0
25	0.0
26	0.0
27	0.0
28	0.0
29	0.0
30	0.0
31	0.0
32	0.0
33	0.0
34	0.0
35	0.0
36	0.0
37	0.0
38	0.0
39	0.3333333333333333
40	0.3333333333333333
41	0.0
42	0.0
43	0.0
44	0.0
45	0.0
46	0.0
47	0.3333333333333333
48	0.3333333333333333
49	0.6666666666666666
50	0.0
51	1.0
52	0.3333333333333333
53	0.0
54	0.0
55	0.0
56	0.0
57	0.0
58	0.3333333333333333
59	0.6666666666666666
60	0.0
61	0.0
62	0.3333333333333333
63	0.0
64	0.0
65	0.0
66	0.0
67	0.0
68	0.0
69	0.0
70	0.0
71	0.0
72	0.0
73	0.0
74	0.0
75	0.0
76	0.0
77	0.0
78	0.0
79	0.0
80	0.0
81	0.0
82	0.0
83	0.0
84	0.0
85	0.0
86	0.0
87	0.0
88	0.0
89	0.0
90	0.0
91	0.0
92	0.0
93	0.0
94	0.0
95	0.0
96	0.0
97	0.0
98	0.0
99	0.0
100	0.0
>>END_MODULE
>>Per base N content	warn
#Base	N-Count
1	0.0
2	0.0
3	0.0
4	0.0
5	0.0
6	0.0
7	0.0
8	0.0
9	0.0
10-14	0.0
15-19	0.0
20-24	8.0
25-29	0.0
30-34	0.0
35-39	0.0
40-44	0.0
45-49	0.0
50-54	0.0
55-59	0.0
60-64	0.0
65-69	0.0
70-74	0.0
75-79	0.0
80-84	0.0
85-89	0.0
90-94	0.0
95-99	0.0
100-104	0.0
105-109	0.0
110-114	0.0
115-119	0.0
120-124	0.0
125-129	0.0
130-134	0.0
135-139	0.0
140-144	0.0
145-149	0.0
150-154	0.0
155-159	0.0
160-164	0.0
165-169	0.0
170-174	0.0
175-179	0.0
180-184	0.0
185-189	0.0
190-194	0.0
195-199	0.0
200-204	0.0
205-209	0.0
210-214	0.0
215-219	0.0
220-224	0.0
225-229	0.0
230-234	0.0
235-239	0.0
240-244	0.0
245-249	0.0
250	0.0
>>END_MODULE
>>Sequence Length Distribution	pass
#Length	Count
250	10.0
>>END_MODULE
>>Sequence Duplication Levels	pass
#Total Deduplicated Percentage	100.0
#Duplication Level	Percentage of deduplicated	Percentage of total
1	100.0	100.0
2	0.0	0.0
3	0.0	0.0
4	0.0	0.0
5	0.0	0.0
6	0.0	0.0
7	0.0	0.0
8	0.0	0.0
9	0.0	0.0
>10	0.0	0.0
>50	0.0	0.0
>100	0.0	0.0
>500	0.0	0.0
>1k	0.0	0.0
>5k	0.0	0.0
>10k+	0.0	0.0
>>END_MODULE
>>Overrepresented sequences	fail
#Sequence	Count	Percentage	Possible Source
GGTATCCCCCGGCAGTGAGGATGGAGCCATGGTCTGCATCATACTCACCG	1	10.0	No Hit
GAGCAGTCGGGCTCAGCGCTNTGCAAATTCTAGTTAGAAACTCACAGTTC	1	10.0	No Hit
AAAATAATCCTAAAAAATAACCTCTATGCCGCCGAACGCTCCGCCTCTAT	1	10.0	No Hit
GTAGTATTCTCATATCACAAGTCCCCAAACTGCATAAGGTGTGGAGTGGA	1	10.0	No Hit
ATATATTCATCCGCCATTATNAGAGTCCGATTACTTTAGAACAGTGCCGC	1	10.0	No Hit
CATCACTACCGCTCAGGAATNTGACGGCAGTCTTAGCGGCGCTCTAGTGC	1	10.0	No Hit
AGATAGCCGAAGATAAAGAGNTCATAACCGTAAAGGCCAGAGACGAGAAC	1	10.0	No Hit
GTGCAGGGGGTGTGGTCAATCCACACTGTTGCTGAGGTGATTGGGTCTCC	1	10.0	No Hit
TCTCAGATGAGCATGCAGCAGCCCAGACTCGCCCCACGCAGTTTGCCAAC	1	10.0	No Hit
CTTCCCCACGGCCCAGACACAAGAGACGACCTCCATAAATCTTTTAGAGG	1	10.0	No Hit
>>END_MODULE
>>Adapter Content	fail
#Position	Illumina Universal Adapter	Illumina Small RNA Adapter	Nextera Transposase Sequence	SOLID Small RNA Adapter
1	0.0	0.0	0.0	0.0
2	0.0	0.0	0.0	0.0
3	0.0	0.0	0.0	0.0
4	0.0	0.0	0.0	0.0
5	0.0	0.0	0.0	0.0
6	0.0	0.0	0.0	0.0
7	0.0	0.0	0.0	0.0
8	0.0	0.0	0.0	0.0
9	0.0	0.0	0.0	0.0
10-14	0.0	0.0	0.0	0.0
15-19	0.0	0.0	0.0	0.0
20-24	0.0	0.0	0.0	0.0
25-29	0.0	0.0	0.0	0.0
30-34	0.0	0.0	0.0	0.0
35-39	0.0	0.0	0.0	0.0
40-44	0.0	0.0	0.0	0.0
45-49	0.0	0.0	0.0	0.0
50-54	0.0	0.0	0.0	0.0
55-59	0.0	0.0	0.0	0.0
60-64	0.0	0.0	0.0	0.0
65-69	0.0	0.0	0.0	0.0
70-74	0.0	0.0	0.0	0.0
75-79	0.0	0.0	0.0	0.0
80-84	0.0	0.0	0.0	0.0
85-89	0.0	0.0	0.0	0.0
90-94	0.0	0.0	0.0	0.0
95-99	0.0	0.0	0.0	0.0
100-104	0.0	0.0	0.0	0.0
105-109	0.0	0.0	2.0	0.0
110-114	0.0	0.0	10.0	0.0
115-119	0.0	0.0	10.0	0.0
120-124	0.0	0.0	10.0	0.0
125-129	0.0	0.0	10.0	0.0
130-134	0.0	0.0	10.0	0.0
135-139	0.0	0.0	10.0	0.0
140-144	0.0	0.0	16.0	0.0
145-149	0.0	0.0	20.0	0.0
150-154	0.0	0.0	20.0	0.0
155-159	0.0	0.0	20.0	0.0
160-164	0.0	0.0	20.0	0.0
165-169	0.0	0.0	20.0	0.0
170-174	0.0	0.0	20.0	0.0
175-179	0.0	0.0	20.0	0.0
180-184	0.0	0.0	20.0	0.0
185-189	0.0	0.0	20.0	0.0
190-194	0.0	0.0	20.0	0.0
195-199	0.0	0.0	20.0	0.0
200-204	0.0	0.0	20.0	0.0
205-209	0.0	0.0	20.0	0.0
210-214	0.0	0.0	22.0	0.0
215-219	0.0	0.0	30.0	0.0
220-224	0.0	0.0	30.0	0.0
225-229	0.0	0.0	30.0	0.0
230-234	0.0	0.0	30.0	0.0
235-238	0.0	0.0	30.0	0.0
>>END_MODULE
>>Kmer Content	pass
>>END_MODULE
##FastQC	0.11.3
>>Basic Statistics	pass
#Measure	Value
Filename	Illumina_SG_R1.fastq
File type	Conventional base calls
Encoding	Sanger / Illumina 1.9
Total Sequences	10
Sequences flagged as poor quality	0
Sequence length	250
%%GC	51
>>END_MODULE
>>Per base sequence quality	fail
#Base	Mean	Median	Lower Quartile	Upper Quartile	10th Percentile	90th Percentile
1	27.6	0.0	0.0	0.0	0.0	0.0
2	25.6	0.0	0.0	0.0	0.0	0.0
3	22.0	0.0	0.0	0.0	0.0	0.0
4	26.8	0.0	0.0	0.0	0.0	0.0
5	28.4	0.0	0.0	0.0	0.0	0.0
6	31.3	0.0	0.0	0.0	0.0	0.0
7	29.4	0.0	0.0	0.0	0.0	0.0
8	29.7	0.0	0.0	0.0	0.0	0.0
9	31.9	0.0	0.0	0.0	0.0	0.0
10-14	31.639999999999997	0.0	0.0	0.0	0.0	0.0
15-19	34.48	0.0	0.0	0.0	0.0	0.0
20-24	30.339999999999996	0.0	0.0	0.0	0.0	0.0
25-29	34.019999999999996	0.0	0.0	0.0	0.0	0.0
30-34	34.839999999999996	0.0	0.0	0.0	0.0	0.0
35-39	35.6	0.0	0.0	0.0	0.0	0.0
40-44	34.56	0.0	0.0	0.0	0.0	0.0
45-49	36.36	0.0	0.0	0.0	0.0	0.0
50-54	35.279999999999994	0.0	0.0	0.0	0.0	0.0
55-59	35.580000000000005	0.0	0.0	0.0	0.0	0.0
60-64	35.56	0.0	0.0	0.0	0.0	0.0
65-69	34.58	0.0	0.0	0.0	0.0	0.0
70-74	34.67999999999999	0.0	0.0	0.0	0.0	0.0
75-79	34.17999999999999	0.0	0.0	0.0	0.0	0.0
80-84	32.220000000000006	0.0	0.0	0.0	0.0	0.0
85-89	34.22	0.0	0.0	0.0	0.0	0.0
90-94	34.49999999999999	0.0	0.0	0.0	0.0	0.0
95-99	34.3	0.0	0.0	0.0	0.0	0.0
100-104	33.52	0.0	0.0	0.0	0.0	0.0
105-109	32.06	0.0	0.0	0.0	0.0	0.0
110-114	28.2	0.0	0.0	0.0	0.0	0.0
115-119	30.639999999999997	0.0	0.0	0.0	0.0	0.0
120-124	28.840000000000003	0.0	0.0	0.0	0.0	0.0
125-129	28.96	0.0	0.0	0.0	0.0	0.0
130-134	25.860000000000003	0.0	0.0	0.0	0.0	0.0
135-139	27.660000000000004	0.0	0.0	0.0	0.0	0.0
140-144	26.78	0.0	0.0	0.0	0.0	0.0
145-149	23.96	0.0	0.0	0.0	0.0	0.0
150-154	23.9	0.0	0.0	0.0	0.0	0.0
155-159	25.2	0.0	0.0	0.0	0.0	0.0
160-164	22.020000000000003	0.0	0.0	0.0	0.0	0.0
165-169	21.639999999999997	0.0	0.0	0.0	0.0	0.0
170-174	17.72	0.0	0.0	0.0	0.0	0.0
175-179	16.36	0.0	0.0	0.0	0.0	0.0
180-184	15.540000000000001	0.0	0.0	0.0	0.0	0.0
185-189	15.88	0.0	0.0	0.0	0.0	0.0
190-194	19.82	0.0	0.0	0.0	0.0	0.0
195-199	16.419999999999998	0.0	0.0	0.0	0.0	0.0
200-204	17.76	0.0	0.0	0.0	0.0	0.0
205-209	15.319999999999999	0.0	0.0	0.0	0.0	0.0
210-214	16.9	0.0	0.0	0.0	0.0	0.0
215-219	17.080000000000002	0.0	0.0	0.0	0.0	0.0
220-224	14.839999999999998	0.0	0.0	0.0	0.0	0.0
225-229	11.760000000000002	0.0	0.0	0.0	0.0	0.0
230-234	13.74	0.0	0.0	0.0	0.0	0.0
235-239	12.54	0.0	0.0	0.0	0.0	0.0
240-244	8.14	0.0	0.0	0.0	0.0	0.0
245-249	3.7399999999999998	0.0	0.0	0.0	0.0	0.0
250	4.5	0.0	0.0	0.0	0.0	0.0
>>END_MODULE
>>Per tile sequence quality	pass
#Tile	Base	Mean
1101	1	0.0
1101	2	0.0
1101	3	0.0
1101	4	0.0
1101	5	0.0
1101	6	0.0
1101	7	0.0
1101	8	0.0
1101	9	0.0
1101	10-14	0.0
1101	15-19	0.0
1101	20-24	0.0
1101	25-29	0.0
1101	30-34	0.0
1101	35-39	0.0
1101	40-44	0.0
1101	45-49	0.0
1101	50-54	0.0
1101	55-59	0.0
1101	60-64	0.0
1101	65-69	0.0
1101	70-74	0.0
1101	75-79	0.0
1101	80-84	0.0
1101	85-89	0.0
1101	90-94	0.0
1101	95-99	0.0
1101	100-104	0.0
1101	105-109	0.0
1101	110-114	0.0
1101	115-119	0.0
1101	120-124	0.0
1101	125-129	0.0
1101	130-134	0.0
1101	135-139	0.0
1101	140-144	0.0
1101	145-149	0.0
1101	150-154	0.0
1101	155-159	0.0
1101	160-164	0.0
1101	165-169	0.0
1101	170-174	0.0
1101	175-179	0.0
1101	180-184	0.0
1101	185-189	0.0
1101	190-194	0.0
1101	195-199	0.0
1101	200-204	0.0
1101	205-209	0.0
1101	210-214	0.0
1101	215-219	0.0
1101	220-224	0.0
1101	225-229	0.0
1101	230-234	0.0
1101	235-239	0.0
1101	240-244	0.0
1101	245-249	0.0
1101	250	0.0
>>END_MODULE
>>Per sequence quality scores	pass
#Quality	Count
12	1.0
13	0.0
14	0.0
15	0.0
16	1.0
17	0.0
18	1.0
19	0.0
20	0.0
21	0.0
22	1.0
23	0.0
24	0.0
25	0.0
26	0.0
27	1.0
28	0.0
29	1.0
30	1.0
31	1.0
32	0.0
33	2.0
>>END_MODULE
>>Per base sequence content	fail
#Base	G	A	T	C
1	40.0	30.0	10.0	20.0
2	20.0	30.0	40.0	10.0
3	20.0	40.0	40.0	0.0
4	10.0	20.0	20.0	50.0
5	0.0	60.0	30.0	10.0
6	40.0	20.0	10.0	30.0
7	10.0	20.0	40.0	30.0
8	10.0	20.0	30.0	40.0
9	40.0	10.0	0.0	50.0
10-14	26.0	22.0	20.0	32.0
15-19	24.0	46.0	12.0	18.0
20-24	21.73913043478261	26.08695652173913	21.73913043478261	30.434782608695656
25-29	20.0	34.0	18.0	28.000000000000004
30-34	22.0	14.000000000000002	28.000000000000004	36.0
35-39	30.0	34.0	18.0	18.0
40-44	22.0	22.0	30.0	26.0
45-49	28.000000000000004	24.0	24.0	24.0
50-54	34.0	20.0	20.0	26.0
55-59	24.0	24.0	12.0	40.0
60-64	12.0	28.000000000000004	32.0	28.000000000000004
65-69	28.000000000000004	16.0	34.0	22.0
70-74	22.0	26.0	24.0	28.000000000000004
75-79	32.0	16.0	20.0	32.0
80-84	24.0	22.0	22.0	32.0
85-89	14.000000000000002	24.0	34.0	28.000000000000004
90-94	28.000000000000004	28.000000000000004	28.000000000000004	16.0
95-99	28.000000000000004	26.0	10.0	36.0
100-104	20.0	16.0	34.0	30.0
105-109	18.0	18.0	24.0	40.0
110-114	30.0	20.0	18.0	32.0
115-119	18.0	28.000000000000004	30.0	24.0
120-124	26.0	28.000000000000004	18.0	28.000000000000004
125-129	14.000000000000002	26.0	22.0	38.0
130-134	32.0	30.0	22.0	16.0
135-139	22.0	22.0	32.0	24.0
140-144	20.0	36.0	16.0	28.000000000000004
145-149	32.0	18.0	28.000000000000004	22.0
150-154	26.0	26.0	28.000000000000004	20.0
155-159	16.0	22.0	24.0	38.0
160-164	26.0	20.0	36.0	18.0
165-169	30.0	24.0	18.0	28.000000000000004
170-174	28.000000000000004	16.0	24.0	32.0
175-179	24.0	34.0	14.000000000000002	28.000000000000004
180-184	24.0	24.0	20.0	32.0
185-189	22.0	24.0	24.0	30.0
190-194	16.0	14.000000000000002	36.0	34.0
195-199	22.0	30.0	20.0	28.000000000000004
200-204	24.0	30.0	26.0	20.0
205-209	24.0	22.0	26.0	28.000000000000004
210-214	20.0	22.0	26.0	32.0
215-219	18.0	24.0	28.000000000000004	30.0
220-224	18.0	24.0	22.0	36.0
225-229	26.0	24.0	20.0	30.0
230-234	20.0	30.0	18.0	32.0
235-239	28.000000000000004	16.0	18.0	38.0
240-244	28.000000000000004	32.0	22.0	18.0
245-249	22.0	34.0	12.0	32.0
250	30.0	40.0	10.0	20.0
>>END_MODULE
>>Per sequence GC content	fail
#GC Content	Count
0	0.0
1	0.0
2	0.0
3	0.0
4	0.0
5	0.0
6	0.0
7	0.0
8	0.0
9	0.0
10	0.0
11	0.0
12	0.0
13	0.0
14	0.0
15	0.0
16	0.0
17	0.0
18	0.0
19	0.0
20	0.0
21	0.0
22	0.0
23	0.0
24	0.0
25	0.0
26	0.0
27	0.0
28	0.0
29	0.0
30	0.0
31	0.0
32	0.0
33	0.0
34	0.0
35	0.0
36	0.0
37	0.0
38	0.0
39	0.3333333333333333
40	0.3333333333333333
41	0.0
42	0.0
43	0.0
44	0.0
45	0.0
46	0.0
47	0.3333333333333333
48	0.3333333333333333
49	0.6666666666666666
50	0.0
51	1.0
52	0.3333333333333333
53	0.0
54	0.0
55	0.0
56	0.0
57	0.0
58	0.3333333333333333
59	0.6666666666666666
60	0.0
61	0.0
62	0.3333333333333333
63	0.0
64	0.0
65	0.0
66	0.0
67	0.0
68	0.0
69	0.0
70	0.0
71	0.0
72	0.0
73	0.0
74	0.0
75	0.0
76	0.0
77	0.0
78	0.0
79	0.0
80	0.0
81	0.0
82	0.0
83	0.0
84	0.0
85	0.0
86	0.0
87	0.0
88	0.0
89	0.0
90	0.0
91	0.0
92	0.0
93	0.0
94	0.0
95	0.0
96	0.0
97	0.0
98	0.0
99	0.0
100	0.0
>>END_MODULE
>>Per base N content	warn
#Base	N-Count
1	0.0
2	0.0
3	0.0
4	0.0
5	0.0
6	0.0
7	0.0
8	0.0
9	0.0
10-14	0.0
15-19	0.0
20-24	8.0
25-29	0.0
30-34	0.0
35-39	0.0
40-44	0.0
45-49	0.0
50-54	0.0
55-59	0.0
60-64	0.0
65-69	0.0
70-74	0.0
75-79	0.0
80-84	0.0
85-89	0.0
90-94	0.0
95-99	0.0
100-104	0.0
105-109	0.0
110-114	0.0
115-119	0.0
120-124	0.0
125-129	0.0
130-134	0.0
135-139	0.0
140-144	0.0
145-149	0.0
150-154	0.0
155-159	0.0
160-164	0.0
165-169	0.0
170-174	0.0
175-179	0.0
180-184	0.0
185-189	0.0
190-194	0.0
195-199	0.0
200-204	0.0
205-209	0.0
210-214	0.0
215-219	0.0
220-224	0.0
225-229	0.0
230-234	0.0
235-239	0.0
240-244	0.0
245-249	0.0
250	0.0
>>END_MODULE
>>Sequence Length Distribution	pass
#Length	Count
250	10.0
>>END_MODULE
>>Sequence Duplication Levels	pass
#Total Deduplicated Percentage	100.0
#Duplication Level	Percentage of deduplicated	Percentage of total
1	100.0	100.0
2	0.0	0.0
3	0.0	0.0
4	0.0	0.0
5	0.0	0.0
6	0.0	0.0
7	0.0	0.0
8	0.0	0.0
9	0.0	0.0
>10	0.0	0.0
>50	0.0	0.0
>100	0.0	0.0
>500	0.0	0.0
>1k	0.0	0.0
>5k	0.0	0.0
>10k+	0.0	0.0
>>END_MODULE
>>Overrepresented sequences	fail
#Sequence	Count	Percentage	Possible Source
GGTATCCCCCGGCAGTGAGGATGGAGCCATGGTCTGCATCATACTCACCG	1	10.0	No Hit
GAGCAGTCGGGCTCAGCGCTNTGCAAATTCTAGTTAGAAACTCACAGTTC	1	10.0	No Hit
AAAATAATCCTAAAAAATAACCTCTATGCCGCCGAACGCTCCGCCTCTAT	1	10.0	No Hit
GTAGTATTCTCATATCACAAGTCCCCAAACTGCATAAGGTGTGGAGTGGA	1	10.0	No Hit
ATATATTCATCCGCCATTATNAGAGTCCGATTACTTTAGAACAGTGCCGC	1	10.0	No Hit
CATCACTACCGCTCAGGAATNTGACGGCAGTCTTAGCGGCGCTCTAGTGC	1	10.0	No Hit
AGATAGCCGAAGATAAAGAGNTCATAACCGTAAAGGCCAGAGACGAGAAC	1	10.0	No Hit
GTGCAGGGGGTGTGGTCAATCCACACTGTTGCTGAGGTGATTGGGTCTCC	1	10.0	No Hit
TCTCAGATGAGCATGCAGCAGCCCAGACTCGCCCCACGCAGTTTGCCAAC	1	10.0	No Hit
CTTCCCCACGGCCCAGACACAAGAGACGACCTCCATAAATCTTTTAGAGG	1	10.0	No Hit
>>END_MODULE
>>Adapter Content	fail
#Position	Illumina Universal Adapter	Illumina Small RNA Adapter	Nextera Transposase Sequence	SOLID Small RNA Adapter
1	0.0	0.0	0.0	0.0
2	0.0	0.0	0.0	0.0
3	0.0	0.0	0.0	0.0
4	0.0	0.0	0.0	0.0
5	0.0	0.0	0.0	0.0
6	0.0	0.0	0.0	0.0
7	0.0	0.0	0.0	0.0
8	0.0	0.0	0.0	0.0
9	0.0	0.0	0.0	0.0
10-14	0.0	0.0	0.0	0.0
15-19	0.0	0.0	0.0	0.0
20-24	0.0	0.0	0.0	0.0
25-29	0.0	0.0	0.0	0.0
30-34	0.0	0.0	0.0	0.0
35-39	0.0	0.0	0.0	0.0
40-44	0.0	0.0	0.0	0.0
45-49	0.0	0.0	0.0	0.0
50-54	0.0	0.0	0.0	0.0
55-59	0.0	0.0	0.0	0.0
60-64	0.0	0.0	0.0	0.0
65-69	0.0	0.0	0.0	0.0
70-74	0.0	0.0	0.0	0.0
75-79	0.0	0.0	0.0	0.0
80-84	0.0	0.0	0.0	0.0
85-89	0.0	0.0	0.0	0.0
90-94	0.0	0.0	0.0	0.0
95-99	0.0	0.0	0.0	0.0
100-104	0.0	0.0	0.0	0.0
105-109	0.0	0.0	2.0	0.0
110-114	0.0	0.0	10.0	0.0
115-119	0.0	0.0	10.0	0.0
120-124	0.0	0.0	10.0	0.0
125-129	0.0	0.0	10.0	0.0
130-134	0.0	0.0	10.0	0.0
135-139	0.0	0.0	10.0	0.0
140-144	0.0	0.0	16.0	0.0
145-149	0.0	0.0	20.0	0.0
150-154	0.0	0.0	20.0	0.0
155-159	0.0	0.0	20.0	0.0
160-164	0.0	0.0	20.0	0.0
165-169	0.0	0.0	20.0	0.0
170-174	0.0	0.0	20.0	0.0
175-179	0.0	0.0	20.0	0.0
180-184	0.0	0.0	20.0	0.0
185-189	0.0	0.0	20.0	0.0
190-194	0.0	0.0	20.0	0.0
195-199	0.0	0.0	20.0	0.0
200-204	0.0	0.0	20.0	0.0
205-209	0.0	0.0	20.0	0.0
210-214	0.0	0.0	22.0	0.0
215-219	0.0	0.0	30.0	0.0
220-224	0.0	0.0	30.0	0.0
225-229	0.0	0.0	30.0	0.0
230-234	0.0	0.0	30.0	0.0
235-238	0.0	0.0	30.0	0.0
>>END_MODULE
>>Kmer Content	pass
>>END_MODULE
""",
}

FASTQ_SCREEN_V0_9_2 = {
    'screen.txt': """#Fastq_screen version: 0.9.2	#Aligner: bowtie	#Reads in subset: 1000000
Genome	#Reads_processed	#Unmapped	%Unmapped	#One_hit_one_genome	%One_hit_one_genome	#Multiple_hits_one_genome	%Multiple_hits_one_genome	#One_hit_multiple_genomes	%One_hit_multiple_genomes	Multiple_hits_multiple_genomes	%Multiple_hits_multiple_genomes
hg19	203860	202199	99.18	2	0.00	0	0.00	136	0.07	1523	0.75
mm9	203860	201733	98.96	89	0.04	62	0.03	58	0.03	1918	0.94
rn4	203860	201822	99.00	48	0.02	21	0.01	93	0.05	1876	0.92
dm3	203860	203844	99.99	0	0.00	0	0.00	0	0.00	16	0.01
ws200	203860	203853	100.00	0	0.00	0	0.00	1	0.00	6	0.00
ecoli	203860	203860	100.00	0	0.00	0	0.00	0	0.00	0	0.00
saccer	203860	203521	99.83	0	0.00	1	0.00	79	0.04	259	0.13
PhiX	203860	203860	100.00	0	0.00	0	0.00	0	0.00	0	0.00
Vectors	203860	203860	100.00	0	0.00	0	0.00	0	0.00	0	0.00
SpR6	203860	203860	100.00	0	0.00	0	0.00	0	0.00	0	0.00

%Hit_no_genomes: 98.92
"""
}

FASTQ_STRAND_V0_0_4 = { 'fastq_strand.txt': """#fastq_strand version: 0.0.4	#Aligner: STAR	#Reads in subset: 10000
#Genome	1st forward	2nd reverse
human	5.17	103.76
"""
}

# Base64 encoded generic example PNG
BASE64_PNG_DATA = "iVBORw0KGgoAAAANSUhEUgAAAEwAAAAqCAIAAACWWCo7AAABFElEQVR4nO2asW3DMBBFvwwv4YZ1BjDgwn1qNlqDN8vnGirEOjO4j2sVUp0JmIJQQwOOAyOxfbwH4YMQVOjdB0EV6khCPSRzzuM8ak2Sm0cP+T8wSTXUexLQd20r50i63k3DJBJyRlqS3/mSMcYQQlrSNEyud+XOQ4r5Fcej1JIiAQIA45y6zgMewJpBZF1LeTz8/Uteg7xpyrUkGUuTgHv+JkV+nvLpJNsgAhEPAPBAfKlvg1ua3O8zql36DMf3nckLNmTMGfNnKgng8IWXzkuaOCdNUgv1EaKA/t3t3vxyTiWHj6mJJk1SCyaphVXyrDmbalI1JqkFk9SCSWrBJLVgklowSS2YpBZMUgtNSHYt/Kz0Dck+c8rgvriAAAAAAElFTkSuQmCC"
