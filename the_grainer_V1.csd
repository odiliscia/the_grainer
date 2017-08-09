;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
/*the_grainer*/
/*A CSOUND ENVIRONMENT FOR SPATIAL AND SPECTRAL GRANULATION
By Oscar Pablo Di Liscia
odiliscia@unq.edu.ar
Research Program "Sistemas Temporales y Síntesis Espacial en el Arte Sonoro"
Escuela Universitaria de Artes
UNQ
Argentina
*/
/*
GENERAL DESCRIPION:
Presently, the environment includes three Csound instruments:
	1-"the_grainer": is the part of the GS environment that creates a
stream of grains computing all its features and calling with the appropriate parameters
another instrument (the_grain).
	2-"the_grain": synthesizes the grains with the parameters computed by its caller 
instrument (the_grainer). The parameters are invariant over the duration of each grain.
	3-"greverb": produces 3D reverberation for all the streams of grains generated
using fast convolution.

PARAMETERS EXPLANATION:

1-"the_grainer":

1.1 Parameters related to GS Synthesis
	grain gap (secs.)			
	grain gap random deviation (secs.)
	grain duration (secs.)			
	grain duration random deviation (secs.)
	grain amplitude 			
	grain gap random deviation
	grain frequency (1=invariant, 0.5= octave lower, 2=octave higher, etc.)			
	grain frequency random deviation (1=invariant, 0.5= octave lower, 2=octave higher, etc.)
	grain starting read point in audio source file (secs.)			
	grain starting read point in audio source file random deviation (secs.)
	grain source audio file (audio file number, see the macro FILESEL)
	grain envelope function  (Csound Table number)

1.2 Parameters related to Spatialisation (see table of room features for spat3di and change values if needed) 			
	grain azimuth angle (normalized from 0 to 1, equals 0 to twoPI, being 0 at the right of the listener)			
	grain azimuth angle random deviation (normalized from 0 to 1, equals 0 to twoPI, being 0 at the right of the listener)	
	grain distance (mts.)			
	grain distance random deviation (mts)
	grain elevation angle (normalized from 0 to 1, equals 0 to twoPI, being 0 at the middle)			
	grain elevation angle random deviation (normalized from 0 to 1, equals 0 to twoPI, being 0 at the middle)

1.3 Parameters related to Spectral treatment (see FFT parameters and change values if needed) 			
	grain bin offset 		(from 0 to (FFTsize/2)-1)			
	grain bin increment 		(from 0 to (FFTsize/2)-1)
	grain bins to synthesize 	(from 0 to (FFTsize/2)-1)
NB1: the bins to synthesize are dependent on the offset and increment. Th present code do the best to
ensure not to overpass the maximal number of bins allowed in each case.
NB2: for offset=0, increment=1 and nbins=(FFTsize/2)-1 the spectral treatmet will be skipped to speed up the computation.

For each one of the aforementioned parameters, there is the possibility of setting
a base value plus a random deviation value that may or may not be invariant over the grain stream. 
There are, at present, four ways (that are handled by the macro GETVAL) of setting these two values: 
a)(Mode 0):Setting base values and random deviation values that are invariant over the stream of
grains.  
b)(Mode 1): Setting a table which will be read to obtain the base and/or the random
deviation values change over the stream of grains. 
c)(Mode 2): Setting a table whose stored values will be accessed through an index that will be generated randomly
to obtain the base and/or the random deviation values.
d)(Mode 3): Setting a random selection between lower and higher values given. 

Each one of the parameters takes 6 pfields of the score, being the first value the "mode".
Examples:

$GETVAL(0'v2'v3'v4'v5'v6'result'ndx)
Mode 0
Takes only the single value of v3 all other argument values are ignored.

$GETVAL(1'v2'v3'v4'v5'v6'result'ndx)
Mode 1
Gets a value from the table "v2" according to index "ndx" (normalized). The value is obtained by
linear interpolation (using tablei), scaled by "v3", added and offset of "v4".
"v5" is the starting read point in the table (normalized) and "v6" is the number of reads that
the table will fit the stream of grains (i.e., the p3 or note duration of "the_grainer"). Negative
values of "v6" are allowed and will produce backwards read. 

$GETVAL(2'v2'v3'v4'v5'v6'result'ndx)
Mode 2
Gets a value from the table "v2" according to an index generated randomly between "v5" and "v6".
"v3" and "v4" are used respectively to scale and offeset the obtained value as well.

$GETVAL(3'v2'v3'v4'v5'v6'result'ndx)
Mode 3
Gets a value generated randomly between "v5" and "v6" all other argument values are ignored.

For all these four methods of getting values there are corresponding score macros:
#define VAL(base) # 0 0 $base 0 0 0 #
#define TAC(tbl'scal'off'pha'ntims) # 1 $tbl $scal $off $pha $ntims #
#define TRA(tbl'scal'off'from'to) # 2 $tbl $scal $off $from $to #
#define RAN(from'to) # 3 0 0 0 $from $to # 
	
*/
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
<CsOptions>

</CsOptions>
<CsInstruments>

sr 	= 44100
ksmps 	= 32	
nchnls 	= 2	;change nchnls to get different output types
		;2=UHJ stereo
		;4=First Order Ambisonic B-format
		;9=Second Order Ambisonic B-format
		;16=Second Order Ambisonic B-format
0dbfs 	= 1

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;GLOBAL VARIABLES
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
giamp		init 	1		;just a master volume control
gipi 		init 	4.*taninv(1.)	;the honorable PI and his father
gipi2		init	2.*gipi
gimaxdel	init	.5		;max delay allowed
/*global audio array for the four output types*/
gaBform[] 	init 	(nchnls <=4? 4 : nchnls) ;up to 3rd order Ambisonics
garvsend	init	0		;global reverb input
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;MACROS
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
/*output types*/
#define STE	#2#
#define AMB1	#4#
#define AMB2	#9#
#define AMB3	#16#
;FFT parameters
#define FFTSIZE  #1024#
#define FFTOVRL  #256#
#define FFTWSIZ  #2048#
#define FFTWTYP  #0#
/*impulse response files, change accordingly, must be a B-Format First order Ambisonic IR*/
#define RVW	#"./IR/York_Centre_W_50ms.wav"#
#define RVX	#"./IR/York_Centre_X_50ms.wav"#
#define RVY	#"./IR/York_Centre_Y_50ms.wav"#
#define RVZ	#"./IR/York_Centre_Z_50ms.wav"#
#define IRDUR	#1.5# /*duration of the impulse response, change accordingly*/
;TABLES WITH ROOM PARAMETERS FOR SPAT3D ARE ALWAYS 99 AND 100
gispatparam	init 0
#define	WL	#	1	#
#define	CFL	#	.5	#
#define	FC	#	3000	#
#define	BP	#	0	#
#define	LP	#	2	#
#define	HP	#	1	#
#define RNDV	#	.1	#
#define	NOE	#0#
#define	FO	#1#
#define	SO	#2#
#define FLTL	#.5#
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
/*for rectangular to polar conversion*/ 
#define POL2REC(az'di'el)		
#
idis		=$di
ix		=cos($az)*cos($el)*$di	
iy		=sin($az)*cos($el)*$di
iz 		=sin($el)*$di
#
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
/*conversion to degrees in Ambisonics coordinates (for the azimuth only)*/
#define AMB_ANGLE(angle'az)		
#					
	$angle = ($az - 0.25) * 360.
	if($angle < 0.) then
		$angle += 360
	endif
#
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
/*UHJ encoding*/
#define UHJ(W'X'Y)
#
aWre, aWim      hilbert $W
aXre, aXim      hilbert $X
aYre, aYim      hilbert $Y
aWXr    =  0.0928*aXre + 0.4699*aWre
aWXiYr  =  0.2550*aXim - 0.1710*aWim + 0.3277*aYre
aleft   =  aWXr + aWXiYr
aright  =  aWXr - aWXiYr
#
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
/*clean the content of an audio rate array*/
;clean (initialize) audio array of 16 values
#define CLEAN_ARRAY16(varr'label)
#
	$varr[0]=0
	$varr[1]=0
	$varr[2]=0
	$varr[3]=0
if nchnls <= $AMB1 goto $label
	$varr[4]=0
	$varr[5]=0
	$varr[6]=0
	$varr[7]=0
	$varr[8]=0
if nchnls == $AMB2 goto $label
	$varr[9]=0
	$varr[10]=0
	$varr[11]=0
	$varr[12]=0
	$varr[13]=0
	$varr[14]=0
	$varr[15]=0
$label:
#
/*replace the first four values of an audio array by other four values*/
#define ADD2ARR4V(varray'v1'v2'v3'v4)
#
$varray[0] = $v1
$varray[1] = $v2
$varray[2] = $v3
$varray[3] = $v4
#
/*add four values to the first four values of an audio array*/
#define ACC2ARR4V(varray'v1'v2'v3'v4)
#
$varray[0] = $varray[0] + $v1
$varray[1] = $varray[1] + $v2
$varray[2] = $varray[2] + $v3
$varray[3] = $varray[3] + $v4
#
/*accumulate up to 16 values of an audio array in other*/
#define ADD2ARR16V(varr1'varr2'label)
#
$varr1[0] = $varr1[0]  + $varr2[0] 
$varr1[1] = $varr1[1]  + $varr2[1]
$varr1[2] = $varr1[2]  + $varr2[2]
$varr1[3] = $varr1[3]  + $varr2[3]
if nchnls <= $AMB1 goto $label
$varr1[4] = $varr1[4]  + $varr2[4]
$varr1[5] = $varr1[5]  + $varr2[5]
$varr1[6] = $varr1[6]  + $varr2[6]
$varr1[7] = $varr1[7]  + $varr2[7]
$varr1[8] = $varr1[8]  + $varr2[8]
if nchnls == $AMB2 goto $label
$varr1[9] = $varr1[9]  + $varr2[9]
$varr1[10] = $varr1[10]  + $varr2[10]
$varr1[11] = $varr1[11]  + $varr2[11]
$varr1[12] = $varr1[12]  + $varr2[12]
$varr1[13] = $varr1[13]  + $varr2[13]
$varr1[14] = $varr1[14]  + $varr2[14]
$varr1[15] = $varr1[15]  + $varr2[15]
$label:
#
/*do output according nchnls*/
#define DO_OUTPUT(outarr'scal'label_1'label_2)
#
if ($scal !=0) goto $label_1
	if (nchnls==$STE) then
	$UHJ($outarr[0]'$outarr[1]'$outarr[2])
		outs aleft, aright
	elseif (nchnls==$AMB1) then 	
		outq	$outarr[0] ,$outarr[1] ,$outarr[2] ,$outarr[3] 
	elseif (nchnls==$AMB2) then 	
		outch	1,$outarr[0] ,2,$outarr[1] ,3,$outarr[2] ,4,$outarr[3] ,5,$outarr[4] ,6,$outarr[5] ,7,$outarr[6] ,8,$outarr[7] ,9,$outarr[8]
	elseif (nchnls==$AMB3) then 	
		outch	1,$outarr[0]  ,2,$outarr[1]   ,3,$outarr[2]   ,4,$outarr[3]   ,5,$outarr[4]   ,6,$outarr[5]   ,7,$outarr[6]  ,8,$outarr[7] ,9,$outarr[8],
		10,$outarr[9] ,11,$outarr[10] ,12,$outarr[11] ,13,$outarr[12] ,14,$outarr[13] ,15,$outarr[14] ,16,$outarr[15]
endif
goto $label_2

$label_1:
	if (nchnls==$STE) then
	$UHJ($outarr[0]'$outarr[1]'$outarr[2])
		outs aleft*$scal, aright*$scal
	elseif (nchnls==$AMB1) then 	
		outq	$outarr[0]*$scal ,$outarr[1]*$scal ,$outarr[2]*$scal ,$outarr[3]*$scal 
	elseif (nchnls==$AMB2) then 	
		outch	1,$outarr[0]*$scal ,2,$outarr[1]*$scal ,3,$outarr[2]*$scal ,4,$outarr[3]*$scal ,5,$outarr[4]*$scal ,6,$outarr[5]*$scal ,7,$outarr[6]*$scal ,8,$outarr[7]*$scal ,9,$outarr[8]*$scal
	elseif (nchnls==$AMB3) then 	
		outch	1,$outarr[0]*$scal  ,2,$outarr[1]*$scal   ,3,$outarr[2]*$scal   ,4,$outarr[3]*$scal   ,5,$outarr[4]*$scal   ,6,$outarr[5]*$scal   ,7,$outarr[6]*$scal  ,8,$outarr[7]*$scal ,9,$outarr[8]*$scal,
		10,$outarr[9]*$scal ,11,$outarr[10]*$scal ,12,$outarr[11]*$scal ,13,$outarr[12]*$scal ,14,$outarr[13]*$scal ,15,$outarr[14]*$scal ,16,$outarr[15]*$scal
endif
$label_2:
#
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
/*select source audio files based on a numeric argument*/
#define FILESEL(arg)
#
if $arg == 1 then
	Sfilename = "./waves/BCL-A3.wav"	;bass clarinet
elseif $arg == 2 then 
	Sfilename = "./waves/BCL-Cs2.wav"	;bass clarinet
elseif $arg == 3 then 
	Sfilename = "./waves/NO_f1_sil3.wav"	;female speech
elseif $arg == 4 then 
	Sfilename = "./waves/NO_f2_gri.wav"	;female shouting
endif
# 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
/*this macro handle the different way of getting grains parameters*/
#define GETVAL(v1'v2'v3'v4'v5'v6'result'ndx)
#
if($v1==0) then 			;single value mode
	$result=$v3
/*v1=mode v2=table v3=scaler v4=offset v5=phase v6=ntimes*/ 
elseif($v1==1) then			;continous table scan
	iaux= $ndx*$v6
	while (iaux < 0) do ;needed for negative indexes
		iaux +=1
	od
	$result table3 iaux, $v2, 1, $v5, 1
	$result*=$v3
	$result+=$v4
/*v1=mode v2=table v3=scaler v4=offset v5=from v6=to*/ 
elseif($v1==2) then			;random table scan
	iaux random $v5, $v6
	$result table abs(int(iaux)), $v2, 0, 0, 1
	$result*=$v3
	$result+=$v4
/*v1=mode v2=na v3=na v4=na v5=from v6=to*/ 
elseif($v1==3) then			;random selection
	$result random $v5, $v6
endif
#
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;TABLE GENERATION FOR spat3di, change accordingly your preferences.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
gispat    	ftgen 	0, 0, 64, -2, $SO, 3, -1, -1, -1, -1,	
/*6 */1,	2,	$RNDV,	$CFL,	$FC,	$FLTL,	.7071,	$LP,	/*ceil*/	
/*14*/1,	2,	$RNDV,	$CFL,	$FC,	$FLTL,	.7071,	$LP,	/*floor*/
/*22*/1,	15,	$RNDV,	$WL,	$FC,	$FLTL,	.7071,	$LP,	/*front*/	
/*30*/1,	15,	$RNDV,	$WL,	$FC,	$FLTL,	.7071,	$LP,	/*back*/
/*38*/1,	15,	$RNDV,	$WL,	$FC,	$FLTL,	.7071,	$LP,	/*right*/
/*46*/1,	15,	$RNDV,	$WL,	$FC,	$FLTL,	.7071,	$LP	/*left*/
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;TABLE GENERATION FOR THE 25 CRITICAL BAND EDGES
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
gicbe    	ftgen 	0, 0, 32, -2,    0.0, 100.0, 200.0, 300.0, 400.0, 510.0, 630.0, 770.0, 920.0, 1080.0,1270.0, 1480.0, 1720.0, 
		2000.0, 2320.0, 2700.0, 3150.0, 3700.0, 4400.0, 5300.0, 6400.0, 7700.0, 9500.0, 12000.0, 15500.0, 20000.0 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;INSTRUMENTS DEFINITIONS
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
/*********************************************************************************************************************/
/*"the_grainer" generates sequences of spatialised and -posssibly- spectrally modified grains
It does so by calling recursively the next instrument ("the_grain") with -possibly- different
parameters at each call.
*/
/*********************************************************************************************************************/
instr the_grainer

indx	=0	;index for table reading	
igacum	=0	;auxiliar variable to store the gap times accumulation
seed	p130	;seed the random number generator

;RECURSIVELY COMPUTE THE PARAMETERS FOR THE GRAINS AND SYNTHESIZE THEM BY CALLING the_grain INSTRUMENT
while (indx < 1) do
	;grains DUR computing
	$GETVAL(p16'p17'p18'p19'p20'p21'idurbas'indx)
	$GETVAL(p22'p23'p24'p25'p26'p27'idurran'indx)
	ibirnd	random -idurran, idurran
	idur	= abs(idurbas + ibirnd)	
	;grains AMP computing
	$GETVAL(p28'p29'p30'p31'p32'p33'iampbas'indx)
	$GETVAL(p34'p35'p36'p37'p38'p39'iampran'indx)
	ibirnd	random -iampran, iampran
	iamp	= abs(iampbas + ibirnd)	
	;grains FREQ computing
	$GETVAL(p40'p41'p42'p43'p44'p45'ifrebas'indx)
	$GETVAL(p46'p47'p48'p49'p50'p51'ifreran'indx)
	ibirnd	random -ifreran, ifreran
	ifreq	= ifrebas + ibirnd 
	;grains ACCESS TIMES in the audio table computing
	$GETVAL(p52'p53'p54'p55'p56'p57'itimbas'indx)
	$GETVAL(p58'p59'p60'p61'p62'p63'itimran'indx)
	ibirnd	random -itimran, itimran
	itim	= abs(itimbas + ibirnd)
	;grains AUDIO TABLE to be used 
	$GETVAL(p64'p65'p66'p67'p68'p69'itau'indx)
	;grains ENVELOPE TABLE to be used
	$GETVAL(p70'p71'p72'p73'p74'p75'iten'indx)
	;grains AZIMUTH ANGLE computing
	$GETVAL(p76'p77'p78'p79'p80'p81'iazibas'indx)
	$GETVAL(p82'p83'p84'p85'p86'p87'iaziran'indx)
	ibirnd	random -iaziran, iaziran
	iazi	= abs(iazibas + ibirnd)
	;grains DISTANCE computing
	$GETVAL(p88'p89'p90'p91'p92'p93'idisbas'indx)
	$GETVAL(p94'p95'p96'p97'p98'p99'idisran'indx)
	ibirnd	random -idisran, idisran
	idis	= abs(idisbas + ibirnd)
	;grains ELEVATION ANGLE computing
	$GETVAL(p100'p101'p102'p103'p104'p105'ielebas'indx)
	$GETVAL(p106'p107'p108'p109'p110'p111'ieleran'indx)
	ibirnd	random -ieleran, ieleran
	iele	= abs(ielebas + ibirnd)
	;grains PHASE VOCODER PARAMETERS computing
	;bin offset
	$GETVAL(p112'p113'p114'p115'p116'p117'ibinoff'indx)
	ibinoff=int(ibinoff)
	;bin increment	
	$GETVAL(p118'p119'p120'p121'p122'p123'ibininc'indx)
	ibininc=int(ibininc)
	;number of bins
	$GETVAL(p124'p125'p126'p127'p128'p129'ibins'indx)
	ibins=int(ibins)
	
	/*CALL "the_grain" INSTRUMENT (SEE ABOVE) TO SYNTHESIZE THE GRAINS*/
	schedule "the_grain", igacum, idur, iamp, ifreq, itim, itau, iten, iazi, idis, iele, ibinoff, ibininc, ibins
	;next grain GAP computing
	$GETVAL(p4'p5'p6'p7'p8'p9'igapbas'indx)
	$GETVAL(p10'p11'p12'p13'p14'p15'igapran'indx)
	ibirnd	random -igapran, igapran
	igap	= abs(igapbas + ibirnd)
	igacum	+= igap		;accumulate gap times
	indx= igacum / p3	;compute a normalized index (0 to 1) to read the tables, if needed

od

endin
/*********************************************************************************************************************/
/*the_grain generates the spatialised and -possibly- spectrally modified grain
it is meant to be called by the_grainer (the previous instrument), but nothing
prevent the user to call it as a single-note generator*/
/*********************************************************************************************************************/
instr	the_grain

inchnls = (nchnls <= $AMB1? 4 : nchnls)
aBform[] init	inchnls	;audio array for output

iamp	= p4		;peak amplitude
ifreq	= p5 		;frequency warping value (1=invariant)
itim	= p6		;skip time in the audio file
$FILESEL(p7)		;audio file name
iten	= int(p8)	;envelope table number
iazi	= p9		;azimuth (0 to 1 = 0 to 2PI)
$AMB_ANGLE(ialpha'iazi)	;azimuth in degrees
iazi	*= gipi2	;azimuth in radians
idis	= p10		;distance
iele	= p11		;elevation (0 to 1 = 0 to 2PI)
ibeta	= iele*360.	;elevation in degrees
iele	*=gipi2		;elevation in radians
ibinoff	= p12		;first bin
ibininc	= p13		;bin increment
ibins	= p14		;number of bins

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;phase vocoder parameters
;compute maximal number of oscillators allowed
;according fft size, bin offset and bin increment
ifftsize= $FFTSIZE
ioverlap= $FFTOVRL
iwinsize= $FFTWSIZ
iwintype= $FFTWTYP
imaxosc	= (ifftsize / 2)-1
;we try to keep all values in range 
ibinoff = ibinoff > imaxosc || ibinoff < 0? 0 : floor(ibinoff)
ibininc = ibininc > imaxosc || ibininc < 1? 1 : floor(ibininc)
ibins   = ibins   > imaxosc || ibins   < 1? 1 : floor(ibins)
;compute max possible bins to resynthesize
imaxpos	= floor((imaxosc - ibinoff) / ibininc)
;we try to ensure do not trespass the maximal number of bins available
ibins	= ( ibins >= imaxpos ? imaxpos  : ibins)
ibins	= ( ibins <= 0 ? 1  : ibins)
;check if the user request the full audio bands, 
ifull	= (ibinoff==0 && ibins==imaxosc && ibininc==1 ? 1 : 0)

;if spectral processing is needed, then frequency warping
;will be done by pvsadsyn, else it will be done by diskin2
if ifull == 1 then
	ifpv=1
else 
	ifpv =ifreq
	ifreq=1
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;amp gate envelope parameters
iseg	= 0.05*p3
isug	= p3-iseg
ifper	= 1/p3
/*read audio signal from a soundfile with frequency warping, if required*/
aenv	poscil	iamp, ifper, iten, 0 
asig	diskin2	Sfilename, ifreq, itim,0,0,8,0
asig	*=aenv
;if the user requested the full audio bands, 
;then there is no need to use pvsanal & pvsaddsyn and the computation is faster
if ifull == 1 goto nopvoc
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;analize
fsig	pvsanal asig, ifftsize, ioverlap, iwinsize, iwintype
;resynthesize
asig	pvsadsyn fsig, ibins, ifpv, ibinoff, ibininc
;ditto...
nopvoc:
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;spatialize
/*extend p3 according the echoes duration*/
p3	= p3 + gimaxdel
$POL2REC(iazi'idis'iele)	;convert from polar to rectangular for spa3di use
imode		=3		;B format output
iucirc		=1		;distance of the unit circle
agate		linseg iamp,isug,iamp,iseg,0
ain		=asig*agate*giamp
;BF ENCODING (first order) of direct sound & echoes
aW,aX,aY,aZ 	spat3di ain,ix,iy,iz,iucirc,gispat,imode,0
/*the dense reverb input is feed by the direct signal with neither amplitude scaling nor delay*/
garvsend += ain
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
/*Compute second or third order Ambisonics B format only for the direct sound if required*/
if(inchnls <= $AMB1) goto ambi1			;skip 2nd or 3rd order computing if nchnls <=4
	/*prepare the feed for bformenc1 */
	isdis		= sqrt(ix*ix + iy*iy + iz*iz)	;warning: this asume that the mic is placed in the origin
	iaten 		= 1 / (isdis +.1)		;distance attenuation according spat3di documentation
	isdel		= isdis / 340.			;distance delay, according spat3di documentation
	ainbfenc	delay ain*iaten, isdel		;delay and attenuate source input signal			
	aBform 		bformenc1 ainbfenc, ialpha, ibeta ;2nd or 3rd order enconding
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
ambi1:
/*
place the four 1st order B-Format signals of the echoes & direct signals obtained with spat3di
in the four first cells of the output array 
*/
$ADD2ARR4V(aBform'aW'aX'aY'aZ)
/*accumulate the output in the global array for the possible use of the sfwriter instrument*/
$ADD2ARR16V(gaBform'aBform'accumdone)
/*perform output according nchnls*/
$DO_OUTPUT(aBform'0'la1'la2)
endin

/*********************************************************************************************************************/
/*dense reverberation is achieved via fast partitioned convolution*/
/*output is B-Format 1st order, impulse responses to be used must meet the same requirements*/
/*********************************************************************************************************************/
instr greverb	/*warning, should be called before if sfwriter is used*/
	
irvsend = p4
idelay	= p5
p3	+= $IRDUR

if(idelay != 0) then
	arvin delay garvsend, idelay
else	
	arvin=garvsend
endif
	
agrW pconvolve arvin*irvsend, $RVW
agrX pconvolve arvin*irvsend, $RVX
agrY pconvolve arvin*irvsend, $RVY
agrZ pconvolve arvin*irvsend, $RVZ
	
/*accumulate the output in the global array for the possible use of the sfwriter instrument*/
$ACC2ARR4V(gaBform'agrW'agrX'agrY'agrZ)
;clean the global reverb input for the next pass	
garvsend=0
/*perform the output according nchnls*/
/*reverb is just 1st order B-Format*/
if (nchnls==$STE) then
	$UHJ(agrW'agrX'agrY)
	outs aleft, aright
elseif (nchnls >= $AMB1) then 	
	outch 	1,agrW ,2,agrX ,3,agrY ,4,agrZ
endif


endin
/*********************************************************************************************************************/
/*********************************************************************************************************************/

</CsInstruments>

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;SCORE
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
<CsScore>
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;envelopes tables and their macros
#define fix	#10#
#define	p_cce	#11#
#define	p_cve	#12#
#define	p_cnc	#13#
#define	p_cnv	#14#
#define	f_cnv	#15#
#define	f_cnc	#16#
#define	io1	#17#
#define	io2	#18#
#define	io3	#19#
#define	io4	#20#
#define	oi1	#21#
#define	oi2	#22#
#define	oi3	#23#
#define	oi4	#24#
#define	stat1	#25#
#define	stat2	#26#
#define	stat3	#27#
#define	stat4	#28#
#define	stat5	#29#
#define	f_lin	#30#
#define	p_lin	#31#
#define tri	#32#
#define trap	#33#
#define ltran	#34#

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
f10	0	1025	7	1	1024	1
f11	0	4097	16	1	26	3	1	4070	-3	0	
f12	0	4097	16	1	26	-3	1	4070	3	0	
f13	0	4097	16	0	26	3	1	4070	-2.9	0	
f14	0	4097	16	0	26	-3	1	4070	3	0	
f15	0	4097	16	0	4040	-3	1	56	3	0	
f16	0	4097	16	0	4040	3	1	56	-3	0	
f17	0	513	16	0	256	3	1	256	-3	0	
f18	0	513	16	0	256	-3	1	256	3	0
f19	0	513	16	0	256	-3	1	256	-3	0
f20	0	513	16	0	256	3	1	256	3	0
f21	0	513	16	0	16	-3	1	240	3	.1	240	-3	1	16	3	0
f22	0	513	16	0	16	-3	1	240	-3	.1	240	3	1	16	3	0
f23	0	513	16	0	16	-3	1	240	-3	.1	240	-3	1	16	3	0
f24	0	513	16	0	16	-3	1	240	3	.1	240	3	1	16	3	0
f25	0	513	16	0	8	3	1	496	0	1	8	-3	0
f26	0	513	16	0	16	3	1	480	0	1	16	-3	0
f27	0	513	16	0	32	3	1	448	0	1	32	-3	0
f28	0	513	16	0	64	3	1	384	1	1	64	-3	0
f29	0	513	16	0	128	3	1	256	1	1	128	-3	0
f30	0	1025	7	0	1024	1
f31	0	1025	7	1	1024	0
f32	0	1025	7	0	512	1	512 	0
f33	0	1024	7	0	341	1	342 	1	341	0
f34	0	1024	7	0	256	0	256 	1	256	1	256	0

;TABLE WITH THE TIME LOCATION OF THE PHONEMES (corresponding to the file "NO_f1_sil3.wav")
f50 	0 	32 	-2  2.423 7.553 9.265 10.824 12.729 17.831 21.347 22.897 24.711 26.542 28.228 29.878 31.781 35.217 36.921 42.124 43.927 45.513 47.306 48.838 50.367 52.092

;TABLES WITH BIN EVOLUTION
;for the FFT bins
f51 	0 	1024 	7	1	128	1	384	0	128	0	384	1	
f52 	0 	1024 	7	0	146	1	732	1	146	0
f53 	0 	1024 	7	0	204	1	616	1	146	0
f54 	0 	1024 	7	0	341	1	342	1	341	0
;for the azimuth angle
f55 	0 	1024 	7	.375	146	.5	732	.5	146	.25
f56 	0 	1024 	7	.5	204	.75	616	.75	146	.25
f57 	0 	1024 	7	.625	341	1	342	1	341	1.25


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;MACROS TO HANDLE PARAMETER SELECTION
;mode=0 means no table, the value of scaler is taken as base 
;all other parameters are ignored
#define VAL(base) # 0 0 $base 0 0 0 #
;The same, but including both base value and deviation
;all other parameters are ignored
#define VALDEV(base'dev) # 0 0 $base 0 0 0 0 0 $dev 0 0 0 #
;mode=1 means table read continously: table number, scaler, offset, phase and number of reads fitting p3 duration
#define TAC(tbl'scal'off'pha'ntims) # 1 $tbl $scal $off $pha $ntims #
;mode=2 means table with discrete values randomly accessed, the two last parameters are start and end values
;in the table
#define TRA(tbl'scal'off'from'to) # 2 $tbl $scal $off $from $to #
;mode=3 means values randomly selected, the two last parameters are start and end values
;in the table
#define RAN(from'to) # 3 0 0 0 $from $to # 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;some macros to handle the pvoc parameters selection
;no pvoc
#define NOPVOC # 0 0 0 0 0 0  0 0 1 0 0 0  0 0 511 0 0 0 #
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;compute increment in semitones
#define INC(s) #[2^($s/12)]#
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;output types
#define STE	#2#
#define AMB1	#4#
#define AMB2	#9#
#define AMB3	#16#

#define	Cs2	#73.16#
#define NYQ	#22050#
#define F2BIN(freq'fftsize) # [$freq /($NYQ/$fftsize)] # 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
/*
;just a test of a call of the_grain instrument directly from the score
;"the_grain" 		start 	dur 	amp 	freq 	tim 	tau 	ten 	azi 	di 	ele 	binoff 		bininc 	bins 			
i	"the_grain"	0	.6	1	$INC(4)	.7	2	$p_cnc	.25	1	0	0		1 	$F2BIN(2000'512)
i	"the_grain"	.1	.8	.	$INC(4)	.	.	$io1	.125	.	.	$F2BIN(2000'512)	$F2BIN(2000'512)		
i	"the_grain"	0	1	4	$INC(4)	.	.	$f_cnc	.375	.	.	$F2BIN(4000'512)	100				
e
*/
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;just one example of spectral and spatial segregation using several sets of critical bands
;as they become appart in space, the higher bands (2 and 3) are scaled in amplitude to
;become more audible
;1-0 to 2000		center
;2-2000 to 5300		right<->center
;3-5300 to nyquist	left <->center
/*					;gap	gaprnd		dur durrnd  		amp 	 		amprnd 	freq freqrnd
i	"the_grainer"	0	9	$VALDEV(.15'0)		$VALDEV(.3'0.0)   	$TAC($stat1'1'0'0'1) 	$VAL(0)	$VALDEV($INC(6)'0)								
					;tim			timrnd		taud		tenv	
					$TAC($f_lin'1.712'0'0'1) $VAL(0)	$VAL(2)		$VAL($io2)
					;azi azirnd		dist distrnd	ele elernd			
					$VALDEV(.25'0.05)	$VALDEV(1'0)	$VALDEV(0'0) 
					;ibinoff		ibinincr	inbins			rseed
					$VAL(0)			$VAL(1)		$VAL($F2BIN(2000'1024))	0

					;gap	gaprnd		dur durrnd  		amp 	 		amprnd 	freq freqrnd
i	"the_grainer"	0	9	$VALDEV(.15'0)		$VALDEV(.3'0.0)   	$TAC($tri'.5'1'0'3) 	$VAL(0)	$VALDEV($INC(6)'0)					
					;tim			timrnd		taud		tenv	
					$TAC($f_lin'1.712'0'0'1) $VAL(0)	$VAL(2)		$VAL($io2)
					;azi 			 azirnd		dist distrnd		ele elernd			
					$TAC($tri'-.125'.25'0'3) $VAL(0)	$VALDEV(1'0)		$VALDEV(0'0) 
					;ibinoff		ibinincr	inbins							rseed
					$VAL($F2BIN(2000'512))	$VAL(1)		$VAL([$F2BIN(5300'512)-$F2BIN(2000'512)])		0

					;gap	gaprnd		dur durrnd  		amp 	 		amprnd 	freq freqrnd
i	"the_grainer"	0	9	$VALDEV(.15'0)		$VALDEV(.3'0.0) 	$TAC($tri'3'1'0'3) 	$VAL(0)	$VALDEV($INC(6)'0)				
					;tim			timrnd		taud		tenv	
					$TAC($f_lin'1.712'0'0'1) $VAL(0)	$VAL(2)		$VAL($io2)
					;azi 			 azirnd		dist distrnd		ele elernd			
					$TAC($tri'.125'.25'0'3) $VAL(0)	$VALDEV(1'0)		$VALDEV(0'0)  
					;ibinoff		ibinincr	inbins		rseed
					$VAL($F2BIN(5300'512))	$VAL(1)		$VAL(511)	0

;			start	dur	send	delay						
i	"greverb"	0	7.12	.02	0.0	
e
*/
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
/*
;other example of spectral, time and spatial segregation using several sets of FFT bins
;bins 0 to 511 jump by 4 (0,4,8...)	above the head of the listener	<->	right  & middle elev  	
;bins 1 to 511 jump by 4 (1,5,9...)	above the head of the listener	<->	center & middle elev 	
;bins 2 to 511 jump by 4 (2,6,10...)   	above the head of the listener	<->	left   & middle elev  	
;bins 1 to 511 jump by 4 (3,7,11...)	above the head of the listener	<->	back   & middle elev  	
;as the four streams become appart in space, thei lenght is shortened and the grains do not overlap anymore

					;gap	gaprnd	dur 	 durrnd  			amp 	amprnd 	freq freqrnd
i	"the_grainer"	0	12	$VALDEV(.15'0)	$TAC($io1'.3'.05'0.5'4) $VAL(0)  	$VALDEV(1'0) 	$VALDEV(1'0)					
					;tim			timrnd		taud		tenv	
					$TAC($f_lin'.712'0'0'1) $VAL(0)		$VAL(2)		$VAL($tri)
					;azi 	azirnd		dist distrnd	ele elernd			
					$VALDEV(0'0)		$VALDEV(1'0)	$TAC($f_lin'.25'0'.5'4) $VAL(0) 
					;ibinoff		ibinincr	inbins		rseed
					$VAL(0)			$VAL(4)		$VAL(511)	0

					;gap	gaprnd	dur 	 durrnd  			amp 	amprnd 	freq freqrnd
i	"the_grainer"	0	12	$VALDEV(.15'0)	$TAC($io1'.3'.05'0.5'4) $VAL(0)  	$VALDEV(1'0) 	$VALDEV(1'0)					
					;tim			timrnd		taud		tenv	
					$TAC($f_lin'.712'0'0'1) $VAL(0)		$VAL(2)		$VAL($tri)
					;azi 	azirnd		dist distrnd	ele elernd			
					$VALDEV(.25'0)		$VALDEV(1'0)	$TAC($f_lin'.25'0'.5'4) $VAL(0)  
					;ibinoff		ibinincr	inbins		rseed
					$VAL(1)			$VAL(4)		$VAL(511)	0

					;gap	gaprnd	dur 	 durrnd  			amp 	amprnd 	freq freqrnd
i	"the_grainer"	0	12	$VALDEV(.15'0)	$TAC($io1'.3'.05'0.5'4) $VAL(0)  	$VALDEV(1'0) 	$VALDEV(1'0)					
					;tim			timrnd		taud		tenv	
					$TAC($f_lin'.712'0'0'1) $VAL(0)		$VAL(2)		$VAL($tri)
					;azi 	azirnd		dist distrnd	ele elernd			
					$VALDEV(.5'0)		$VALDEV(1'0)	$TAC($f_lin'.25'0'.5'4) $VAL(0)   
					;ibinoff		ibinincr	inbins		rseed
					$VAL(2)			$VAL(4)		$VAL(511)	0

					;gap	gaprnd	dur 	 durrnd  			amp 	amprnd 	freq freqrnd
i	"the_grainer"	0	12	$VALDEV(.15'0)	$TAC($io1'.3'.05'0.5'4) $VAL(0)  	$VALDEV(1'0) 	$VALDEV(1'0)					
					;tim			timrnd		taud		tenv	
					$TAC($f_lin'.712'0'0'1) $VAL(0)		$VAL(2)		$VAL($tri)
					;azi 	azirnd		dist distrnd	ele elernd			
					$VALDEV(.75'0)		$VALDEV(1'0)	$TAC($f_lin'.25'0'.5'4) $VAL(0)  
					;ibinoff		ibinincr	inbins		rseed
					$VAL(3)			$VAL(4)		$VAL(511)	0
;			start	dur	send	delay						
i	"greverb"	0	7.12	.02	0.0	
e
*/
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
/*
;Time, spectral and Spatial Segregation on a vocal sequence of sounds: wovlels and consonants are differently placed in space.
;Consonants (first note) are placed allways in the centre.
;The wovels are divided in two bands the first one form 0 Hz to 2000 Hz (second note), and the second from (third note)
;2000 to the Nyquist freq., and segregated in space differently as well and placed randomly at the sides.
;amplitudes were adjusted to make the effect more apparent. 
					;gap	gaprnd	dur durrnd     amp amprnd    freq freqrnd
i	"the_grainer"	0	15	$VALDEV(.25'.1)	$VALDEV(.07'0) $VALDEV(4'0)  $VALDEV(1'0.025)					
					;tim		 timrnd	  taud		tenv	
					$TRA(50'1'0'0'21) $VAL(0) $VAL(3) $VAL($p_cnv)
					;azi 	azirnd	dist distrnd	ele elernd			
					$VALDEV(.25'0)	$VALDEV(1'0)	$VALDEV(0'0)  
					;ibinoff ibinincr inbins		rseed
					$NOPVOC					1.55

					;gap	gaprnd	dur durrnd     amp amprnd    freq freqrnd
i	"the_grainer"	0	15	$VALDEV(.25'.1)	$VALDEV(.3'0) $VALDEV(1'0)  $VALDEV(1'0.025)					
					;tim		    timrnd	taud	tenv	
					$TRA(50'1'.3'0'21) $VAL(0) 	$VAL(3)	$VAL($stat2)
					;azi 	azirnd		dist distrnd	ele elernd			
					$VALDEV(.125'.125)	$VALDEV(1'0)	$VALDEV(0'0)  
					;ibinoff ibinincr inbins		rseed
					$VAL(0)	 $VAL(1)  $VAL(46)		.


i	"the_grainer"	0	15	$VALDEV(.25'.1)	$VALDEV(.3'0) $VALDEV(2'0)  $VALDEV(1'0.025)					
					;tim		    timrnd	taud	tenv	
					$TRA(50'1'.3'0'21) $VAL(0) 	$VAL(3)	$VAL($stat2)
					;azi 	azirnd		dist distrnd	ele elernd			
					$VALDEV(.375'.125)	$VALDEV(1'0)	$VALDEV(0'0)  
					;ibinoff ibinincr inbins		rseed
					$VAL(47) $VAL(1)  $VAL(511)		.

i	"greverb"	0	15	.01	0.0
e
*/
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
/*
Distance and spectrum segregation
The spectrum of a speech voice shouting a text is divided in five bands
		bin	Hz.
Band 1:		0-50  	(0-2000)
Band 2:		51-75 	(2000-3000)
Band 3:		76-100 	(3000-4000)
Band 4:		101-126 (4000-5400)
Band 5:		127-511 (5400-nyquist)

Each band changes in distance from 1 to 10 mts. according an exponential function
which is read starting in different points for each one. 
*/
/*	
					;gap	gaprnd	dur durrnd     		amp amprnd    freq freqrnd
i	"the_grainer"	0	15	$VALDEV(.25'.1)	$VALDEV(.5'0) 	$VALDEV(1'0)  $VALDEV(1'0.025)					
					;tim	timrnd	  	taud		tenv	
					$VALDEV(10.75'10.75)   $VAL(4) 		$VAL($stat1)
					;azi 	azirnd		dist distrnd			ele elernd			
					$VALDEV(.25'0)		$TAC($io1'5'1'0'5) $VAL(0)	$VALDEV(0'0)  
					;ibinoff ibinincr inbins		rseed
					$VAL(0) $VAL(1)  $VAL(50)		1.55
			
					;gap	gaprnd	dur durrnd     amp amprnd    freq freqrnd
i	"the_grainer"	0	15	$VALDEV(.25'.1)	$VALDEV(.5'0) $VALDEV(2'0)  $VALDEV(1'0.025)					
					;tim	timrnd	  	taud		tenv	
					$VALDEV(10.75'10.75)   $VAL(4) 		$VAL($stat1)
					;azi 	azirnd		dist distrnd			ele elernd			
					$VALDEV(.25'0)		$TAC($io1'9'1'0.6'5) $VAL(0)	$VALDEV(0'0)  
					;ibinoff ibinincr inbins		rseed
					$VAL(51) $VAL(1)  $VAL(24)		1.55
	
					;gap	gaprnd	dur durrnd     amp amprnd    freq freqrnd
i	"the_grainer"	0	15	$VALDEV(.25'.1)	$VALDEV(.5'0) $VALDEV(3'0)  $VALDEV(1'0.025)					
					;tim	timrnd	  	taud		tenv	
					$VALDEV(10.75'10.75)   $VAL(4) 		$VAL($stat1)
					;azi 	azirnd		dist distrnd			ele elernd			
					$VALDEV(.25'0)		$TAC($io1'9'1'0.2'5) $VAL(0)	$VALDEV(0'0)  
					;ibinoff ibinincr inbins		rseed
					$VAL(76) $VAL(1)  $VAL(24)		1.55	

					;gap	gaprnd	dur durrnd     amp amprnd    freq freqrnd
i	"the_grainer"	0	15	$VALDEV(.25'.1)	$VALDEV(.5'0) $VALDEV(4'0)  $VALDEV(1'0.025)					
					;tim	timrnd	  	taud		tenv	
					$VALDEV(10.75'10.75)   $VAL(4) 		$VAL($stat1)
					;azi 	azirnd		dist distrnd			ele elernd			
					$VALDEV(.25'0)		$TAC($io1'9'1'0.4'5) $VAL(0)	$VALDEV(0'0)  
					;ibinoff ibinincr inbins		rseed
					$VAL(101) $VAL(1)  $VAL(25)		1.55	

					;gap	gaprnd	dur durrnd     amp amprnd    freq freqrnd
i	"the_grainer"	0	15	$VALDEV(.25'.1)	$VALDEV(.5'0) $VALDEV(4'0)  $VALDEV(1'0.025)					
					;tim	timrnd	  	taud		tenv	
					$VALDEV(10.75'10.75)   $VAL(4) 		$VAL($stat1)
					;azi 	azirnd		dist distrnd			ele elernd			
					$VALDEV(.25'0)		$TAC($io1'9'1'0.8'5) $VAL(0)	$VALDEV(0'0)  
					;ibinoff ibinincr inbins		rseed
					$VAL(127) $VAL(1)  $VAL(511)		1.55

i	"greverb"	0	15	.01	0.0
*/

/*
Time, spectrum, pitch and spatial segregation
The spectrum of the source signal is divded in two bands (0-3150Hz and 3150Hz to Nyquist freq.)
The more the gap of the two sequences becomes more irregular, they de-sinchronize and also are
slighty modulated in pitch and moves from the centre to both sides. The more regular is the sequence,
they move to the centre, have no pitch modulation and fuse into a unique sound.
*/
/*
					;gap		gaprnd			dur durrnd     	amp amprnd    freq freqrnd
i	"the_grainer"	0	10	$VAL(.25) 	$TAC($ltran'.1'0'0'2)	$VALDEV(.1'.01) $VALDEV(1'0)  $VAL(1) $TAC($ltran'.01'0'0'2) 					
					;tim	timrnd	  	taud		tenv	
					$VALDEV(.4'.01) 	$VAL(1) 	$VAL($io1)
					;azi 	azirnd				dist distrnd	ele elernd			
					$TAC($ltran'.125'.25'0'2) $VAL(0)	$VALDEV(1'0)	$VALDEV(0'0)  
					;ibinoff 	ibinincr 	inbins			rseed
					$VAL(0)		$VAL(1)		$VAL($F2BIN(3150'512))	1.55

					;gap		gaprnd			dur durrnd     	amp amprnd    freq freqrnd
i	"the_grainer"	0	10	$VAL(.25) 	$TAC($ltran'.1'0'0'2)	$VALDEV(.1'.01) $VALDEV(3'0)  $VAL(1) $TAC($ltran'.01'0'0'2) 					
					;tim	timrnd	  	taud		tenv	
					$VALDEV(.4'.01) 	$VAL(1) 	$VAL($io1)
					;azi 	azirnd				dist distrnd	ele elernd			
					$TAC($ltran'-.125'.25'0'2) $VAL(0)	$VALDEV(1'0)	$VALDEV(0'0)  
					;ibinoff 		ibinincr 	inbins		rseed
					$VAL($F2BIN(3150'512))	$VAL(1)		$VAL(511)	1.6
i	"greverb"	0	15	.01	0.0
e
*/
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;base bands				;gap	gaprnd		dur durrnd     			amp amprnd    	freq freqrnd  
i	"the_grainer"	0	16	$VALDEV(.05'0.01)	$TAC(51'.15'0.05'0'1) $VAL(0)   $VALDEV(.75'0)	$VALDEV(1'0)				
					;tim			timrnd	  	taud		tenv	
					$TAC($f_lin'1'0'0'1) 	$VAL(0)   	$VAL(1) 	$VAL($tri)
					;azi 			azirnd	dist distrnd	ele elernd			
					$TAC($f_lin'1'0'.25'1)	$VAL(0)	$VALDEV(1'0)	$VALDEV(0'0)  			
					;ibinoff 	ibinincr 	inbins			rseed
					$VAL(0)		$VAL(1)		$TAC(51'75'50'0'1)	1.55

/*
					;gap	gaprnd		dur durrnd     			amp amprnd    	freq freqrnd
i	"the_grainer"	0	16	$VALDEV(.05'0.01)	$TAC(51'.002'0'0'1) $VAL(0)   	$VALDEV(.75'0)	$VALDEV(1'0)				
					;tim			timrnd	  	taud		tenv	
					$TAC($f_lin'1'0'0'1) 	$VAL(0)   	$VAL(1) 	$VAL($tri)
					;azi 			azirnd		dist distrnd	ele elernd			
					$TAC($f_lin'1'0'.25'1)	$VAL(0)		$VALDEV(1'0)	$VALDEV(0'0)  			
					;ibinoff 	ibinincr 	inbins		rseed
					$VAL(126)	$VAL(1)		$VAL(384)	1.55
*/
;evolutive bands
					;gap	gaprnd		dur durrnd     			amp amprnd    freq freqrnd
i	"the_grainer"	2	14	$VALDEV(.05'0.01)	$TAC(52'-.15'.2'0'1) $VAL(0) 	$VALDEV(2'0)  $VALDEV(1'0)					
					;tim			timrnd	  	taud		tenv	
					$TAC($f_lin'1'0'0'1) 	$VAL(0)   	$VAL(1) 	$VAL($tri)
					;azi 			azirnd			dist distrnd	ele elernd			
					$TAC(55'1'0'0'1)	$TAC(53'.125'0'0'1) 	$VALDEV(1'0)	$VALDEV(0'0)
					;ibinoff 		ibinincr 	inbins			rseed
					$TAC(52'-25'125'0'1)	$VAL(1)		$TAC(52'24'1'0'1)	1.55

					;gap	gaprnd		dur durrnd     			amp amprnd    freq freqrnd
i	"the_grainer"	4	10	$VALDEV(.05'0.01)	$TAC(53'-.15'.2'0'1) $VAL(0) 	$VALDEV(1.5'0)  $VALDEV(1'0)					
					;tim			timrnd	  	taud		tenv	
					$TAC($f_lin'1'0'0'1) 	$VAL(0)   	$VAL(1) 	$VAL($tri)
					;azi 			azirnd			dist distrnd	ele elernd			
					$TAC(56'1'0'0'1)	$TAC(53'.125'0'0'1) 	$VALDEV(1'0)	$VALDEV(0'0)
					;ibinoff 		ibinincr 	inbins			rseed
					$TAC(53'-25'100'0'1)	$VAL(1)		$TAC(53'24'1'0'1)	1.55

					;gap	gaprnd		dur durrnd     			amp amprnd    freq freqrnd
i	"the_grainer"	6	6	$VALDEV(.05'0.01)	$TAC(54'-.15'.2'0'1) $VAL(0) 	$VALDEV(.75'0)  $VALDEV(1'0)				
					;tim			timrnd	  	taud		tenv	
					$TAC($f_lin'1'0'0'1) 	$VAL(0)   	$VAL(1) 	$VAL($tri)
					;azi 			azirnd			dist distrnd	ele elernd			
					$TAC(57'1'0'0'1)	$TAC(53'.125'0'0'1) 	$VALDEV(1'0)	$VALDEV(0'0)
					;ibinoff 		ibinincr 	inbins			rseed
					$TAC(54'-25'75'0'1)	$VAL(1)		$TAC(54'24'1'0'1)	1.55


</CsScore>


</CsoundSynthesizer>
