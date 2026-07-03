# Encuentros notables en el catálogo Gaia DR3

> **Estado:** 🟡 EN CURSO (frente P1 de `planning/PUBLISH_PUSH_PLAN.md`)
> **Fecha:** 2026-07-02
> **Fuente:** `data/output/encounters_characterized_full.parquet` — 72,236,904 encuentros 3D reales (una fila por par, ya deduplicado por máxima aproximación).
>
> Este catálogo registra la **mínima distancia física en 3D** entre pares de asteroides durante la ventana Gaia DR3 (jul 2014 – may 2017), no co-localizaciones aparentes en el plano del cielo.

## Metodología y limitaciones

- Los diámetros derivan de `H` con albedo por clase cuando no hay medida directa; deben leerse como estimaciones de orden de magnitud, no como diámetros medidos. Los cortes por tamaño son por tanto aproximados.
- La propagación de base del catálogo es Kepler de dos cuerpos; el refinamiento N-cuerpos se aplica sólo al subset de determinación de masas. Las distancias mínimas de esta minería pueden tener sesgo cerca de resonancias o encuentros planetarios. Cualquier evento seleccionado como candidato requiere revalidación N-cuerpos antes de publicar.
- El presupuesto de completitud del catálogo (censura ~0.70 %, recall prefiltro ~76 % en el tail adverso) implica que faltan algunos encuentros genuinos; las tablas de abajo son un piso, no un censo exhaustivo.
- El proxy de familia (§4) es proximidad en elementos osculantes `(a, e, i)`, **no** clasificación en elementos propios. Señala candidatos, no confirma pertenencia a familia.

## 1. Encuentros grande-grande

Ambos cuerpos por encima del umbral de diámetro. Son los más raros: dos cuerpos masivos que se aproximan en 3D. Ordenados por distancia mínima.

### 1.1 Ambos D ≳ 50 km (30 en el top 30)

| number_1 | designation_1 | number_2 | designation_2 | date_utc | dist_km | rel_vel_km_s | diameter_1_km | diameter_2_km | class_1 | class_2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 305 | (305) Gordonia | 830 | (830) Petropolitana | 2014-09-18 | 9.876e+05 | 3.280 | 62.584 | 53.760 | Other | Other |
| 426 | (426) Hippo | 788 | (788) Hohensteina | 2015-04-17 | 1.147e+06 | 8.824 | 73.529 | 64.634 | MBA | MBA |
| 145 | (145) Adeona | 780 | (780) Armenia | 2016-11-09 | 1.505e+06 | 6.339 | 84.035 | 56.294 | MBA | MBA |
| 739 | (739) Mandeville | 415 | (415) Palatia | 2017-04-05 | 1.699e+06 | 8.363 | 70.870 | 51.105 | Other | Other |
| 675 | (675) Ludmilla | 80 | (80) Sappho | 2016-02-15 | 1.710e+06 | 2.969 | 92.995 | 90.045 | MBA | MBA |
| 51 | (51) Nemausa | 91 | (91) Aegina | 2016-07-13 | 2.082e+06 | 4.417 | 120 | 60.598 | Other | Other |
| 46 | (46) Hestia | 482 | (482) Petrina | 2016-07-12 | 2.460e+06 | 4.226 | 75.589 | 60.598 | MBA | Other |
| 758 | (758) Mancunia | 740 | (740) Cantabia | 2016-08-13 | 2.630e+06 | 3.695 | 82.882 | 53.760 | Other | Other |
| 739 | (739) Mandeville | 754 | (754) Malabar | 2016-03-14 | 2.723e+06 | 5.258 | 70.870 | 51.578 | Other | Other |
| 200 | (200) Dynamene | 1304 | (1304) Arosa | 2014-12-25 | 2.846e+06 | 7.877 | 79.152 | 51.341 | Other | Other |
| 81 | (81) Terpsichore | 404 | (404) Arsinoe | 2015-04-15 | 2.854e+06 | 9.008 | 71.526 | 56.035 | Other | MBA |
| 776 | (776) Berbericia | 195 | (195) Eurykleia | 2017-02-09 | 2.875e+06 | 5.973 | 103 | 56.035 | Other | MBA |
| 27 | (27) Euterpe | 116 | (116) Sirona | 2017-02-19 | 2.953e+06 | 1.893 | 141 | 96.931 | Other | MBA |
| 89 | (89) Julia | 335 | (335) Roberta | 2015-01-28 | 2.959e+06 | 7.083 | 170 | 57.340 | Other | Other |
| 275 | (275) Sapientia | 182 | (182) Elsa | 2017-02-10 | 3.185e+06 | 5.019 | 60.320 | 53.267 | MBA | Other |
| 1 | (1) Ceres | 197 | (197) Arete | 2017-05-02 | 3.354e+06 | 4.485 | 763 | 51.816 | MBA | Other |
| 635 | (635) Vundtia | 830 | (830) Petropolitana | 2017-03-30 | 3.418e+06 | 4.124 | 56.035 | 53.760 | MBA | Other |
| 976 | (976) Benjamina | 701 | (701) Oriola | 2014-07-24 | 3.429e+06 | 0.993 | 50.870 | 50.172 | MBA | MBA |
| 23 | (23) Thalia | 481 | (481) Emita | 2014-10-14 | 3.627e+06 | 4.437 | 145 | 61.725 | Other | Other |
| 90 | (90) Antiope | 508 | (508) Princetonia | 2015-11-11 | 3.671e+06 | 4.371 | 78.788 | 74.210 | Other | Other |
| 487 | (487) Venetia | 663 | (663) Gerlinde | 2015-05-20 | 3.874e+06 | 8.131 | 83.649 | 51.105 | Other | MBA |
| 38 | (38) Leda | 796 | (796) Sarita | 2016-02-20 | 3.877e+06 | 9.159 | 76.995 | 53.267 | Other | MBA |
| 53 | (53) Kalypso | 151 | (151) Abundantia | 2015-02-11 | 3.900e+06 | 4.657 | 61.441 | 53.760 | Other | Other |
| 154 | (154) Bertha | 338 | (338) Budrosa | 2015-01-24 | 3.964e+06 | 7.232 | 108 | 70.870 | Other | Other |
| 13 | (13) Egeria | 113 | (113) Amalthea | 2016-02-23 | 4.163e+06 | 5.457 | 159 | 63.454 | Other | MBA |
| 110 | (110) Lydia | 312 | (312) Pierretta | 2016-08-05 | 4.244e+06 | 3.774 | 97.828 | 59.219 | Other | MBA |
| 747 | (747) Winchester | 441 | (441) Bathilde | 2014-09-06 | 4.389e+06 | 8.522 | 103 | 70.544 | Other | MBA |
| 109 | (109) Felicitas | 472 | (472) Roma | 2015-06-11 | 4.461e+06 | 7.756 | 63.163 | 58.406 | MBA | MBA |
| 88 | (88) Thisbe | 774 | (774) Armor | 2017-02-16 | 4.622e+06 | 1.369 | 139 | 61.725 | MBA | MBA |
| 51 | (51) Nemausa | 197 | (197) Arete | 2016-12-14 | 4.686e+06 | 5.971 | 120 | 51.816 | Other | Other |


### 1.2 Ambos D ≳ 100 km (2 en el top 30)

| number_1 | designation_1 | number_2 | designation_2 | date_utc | dist_km | rel_vel_km_s | diameter_1_km | diameter_2_km | class_1 | class_2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | (1) Ceres | 57 | (57) Mnemosyne | 2017-01-10 | 6.515e+06 | 7.251 | 763 | 139 | MBA | MBA |
| 7 | (7) Iris | 44 | (44) Nysa | 2014-08-13 | 7.257e+06 | 5.351 | 281 | 139 | MBA | Other |


## 2. Encuentros extremos

### 2a. Mínima distancia absoluta del catálogo

| number_1 | designation_1 | number_2 | designation_2 | date_utc | dist_km | rel_vel_km_s | diameter_1_km | diameter_2_km | class_1 | class_2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 153222 | (153222) 2000 YD43 | 238587 | (238587) 2004 YX3 | 2016-03-11 | 1094 | 6.898 | 2.821 | 1.780 | Other | Other |
| 15072 | (15072) Landolt | 387599 | (387599) 2001 XF180 | 2014-12-08 | 1760 | 5.867 | 3.239 | 1.624 | Other | Other |
| 270730 | (270730) 2002 QE130 | 366918 | (366918) 2005 UC211 | 2016-05-29 | 2291 | 3.102 | 1.624 | 1.176 | Other | Other |
| 117065 | (117065) 2004 KD9 | 439086 | (439086) 2011 QP5 | 2015-07-27 | 2535 | 8.248 | 3.719 | 1.864 | MBA | Other |
| 52249 | (52249) 1981 EK21 | 408138 | (408138) 2013 CL75 | 2016-11-17 | 2590 | 3.162 | 3.392 | 1.414 | MBA | Other |
| 412792 | (412792) 2014 PU21 | 161150 | (161150) 2002 SL25 | 2016-08-25 | 2936 | 3.985 | 1.350 | 1.073 | Other | Other |
| 17067 | (17067) 1999 GF19 | 236737 | (236737) 2007 JC18 | 2015-02-20 | 2966 | 4.873 | 8.137 | 1.864 | Other | Other |
| 436353 | (436353) 2010 JC112 | 435807 | (435807) 2008 VV60 | 2015-12-19 | 3519 | 6.258 | 1.624 | 1.290 | Other | Other |
| 110273 | (110273) 2001 SX251 | 221390 | (221390) 2005 YS34 | 2015-04-17 | 3802 | 2.549 | 2.821 | 1.864 | MBA | Other |
| 209619 | (209619) 2005 AT19 | 304025 | (304025) 2006 DR59 | 2014-12-22 | 3874 | 4.111 | 1.350 | 1.123 | Other | Other |
| 419812 | (419812) 2010 WL69 | 184372 | (184372) 2005 JH160 | 2015-08-20 | 3894 | 4.348 | 1.123 | 0.852 | Other | Other |
| 108480 | (108480) 2001 KJ60 | 227202 | (227202) 2005 QG107 | 2014-10-18 | 4048 | 4.926 | 1.952 | 1.624 | MBA | Other |
| 247365 | (247365) 2001 XC38 | 173492 | (173492) 2000 SZ187 | 2015-01-28 | 4274 | 5.339 | 2.241 | 1.624 | Other | Other |
| 260384 | (260384) 2004 VP60 | 320613 | (320613) 2008 CW21 | 2015-07-13 | 4379 | 13.886 | 2.573 | 1.952 | Other | Other |
| 97940 | (97940) 2000 QW116 | 225717 | (225717) 2001 RX32 | 2016-12-31 | 4719 | 3.597 | 1.952 | 1.624 | Other | Other |
| 98674 | (98674) 2000 WQ168 | 268752 | (268752) 2006 QX110 | 2016-10-16 | 4727 | 9.120 | 2.573 | 1.864 | Other | Other |
| 190676 | (190676) 2001 BO34 | 355914 | (355914) 2008 XR48 | 2016-02-02 | 4749 | 2.419 | 3.094 | 1.290 | MBA | Other |
| 62903 | (62903) 2000 UK106 | 251913 | (251913) 1999 VJ153 | 2015-01-08 | 4835 | 3.140 | 3.719 | 2.347 | Other | Other |
| 247342 | (247342) 2001 UL210 | 453840 | (453840) 2011 SZ252 | 2017-01-21 | 4963 | 7.278 | 2.044 | 0.892 | Other | Other |
| 2974 | (2974) Holden | 297875 | (297875) 2002 CG123 | 2016-05-20 | 5132 | 5.964 | 7.771 | 1.290 | Other | Other |
| 280677 | (280677) 2005 EL185 | 330942 | (330942) 2009 SL244 | 2017-03-13 | 5156 | 5.202 | 2.954 | 1.232 | Other | Other |
| 168429 | (168429) 1998 SP90 | 199909 | (199909) 2007 GY18 | 2016-01-20 | 5889 | 2.823 | 1.952 | 1.624 | Other | Other |
| 98030 | (98030) 2000 RN7 | 387836 | (387836) 2004 HZ53 | 2016-07-08 | 5916 | 7.146 | 3.895 | 1.864 | MBA | MBA |
| 326673 | (326673) 2002 VB5 | 377477 | (377477) 2005 BV33 | 2017-04-10 | 5922 | 10.054 | 2.241 | 0.777 | Other | Other |
| 173444 | (173444) 2000 LG3 | 120767 | (120767) 1998 BS26 | 2015-01-21 | 5995 | 8.468 | 2.954 | 2.347 | MBA | Other |
| 183009 | (183009) 2002 PZ83 | 215420 | (215420) 2002 GN149 | 2016-12-26 | 6174 | 0.333 | 1.481 | 1.350 | Other | Other |
| 91317 | (91317) 1999 GV17 | 301935 | (301935) 2000 AG66 | 2016-04-24 | 6265 | 1.890 | 2.954 | 2.241 | Other | Other |
| 146405 | (146405) 2001 QJ185 | 295006 | (295006) 2008 EX42 | 2016-01-25 | 6427 | 3.657 | 2.457 | 1.624 | Other | Other |
| 177675 | (177675) 2005 ED98 | 164246 | (164246) 2004 TV113 | 2016-11-28 | 6484 | 5.685 | 1.624 | 1.290 | Other | Other |
| 1823 | (1823) Gliese | 231600 | (231600) 2009 AA19 | 2016-03-26 | 6565 | 5.566 | 10.727 | 1.550 | Other | Other |


### 2b. Encuentros más lentos (candidatos naturales a masa)

Velocidad relativa mínima (con al menos un cuerpo D ≳ 5 km). Una v_rel baja prolonga la interacción gravitatoria ⇒ deflexión mayor y más medible.

| number_1 | designation_1 | number_2 | designation_2 | date_utc | dist_km | rel_vel_km_s | diameter_1_km | diameter_2_km | class_1 | class_2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 3791 | (3791) Marci | 110964 | (110964) 2001 UW170 | 2017-01-18 | 5.847e+06 | 0.006 | 13.504 | 2.457 | Other | Other |
| 47042 | (47042) 1998 WP3 | 99861 | (99861) Tscharnuter | 2014-07-27 | 6.029e+06 | 0.015 | 5.629 | 3.239 | MBA | MBA |
| 7203 | (7203) Sigeki | 73664 | (73664) 1981 EE34 | 2015-02-28 | 7.072e+06 | 0.031 | 8.137 | 1.481 | MBA | Other |
| 3384 | (3384) Daliya | 309262 | (309262) 2007 RZ93 | 2016-12-20 | 4.444e+06 | 0.031 | 6.172 | 1.232 | Other | Other |
| 7479 | (7479) 1994 EC1 | 14954 | (14954) 1996 DL | 2017-05-12 | 7.099e+06 | 0.032 | 6.463 | 3.094 | Other | Other |
| 51906 | (51906) 2001 QN51 | 287055 | (287055) 2002 QR134 | 2015-01-03 | 7.344e+06 | 0.033 | 5.629 | 1.700 | Other | Other |
| 13412 | (13412) Guerrieri | 129102 | (129102) 2004 XO9 | 2015-08-10 | 1.825e+06 | 0.038 | 5.895 | 2.573 | Other | Other |
| 17766 | (17766) 1998 ES3 | 199454 | (199454) 2006 DQ39 | 2016-05-31 | 2.267e+06 | 0.040 | 7.421 | 2.140 | MBA | Other |
| 11907 | (11907) Naranen | 268636 | (268636) 2006 DS101 | 2016-10-13 | 7.080e+06 | 0.041 | 7.771 | 1.864 | MBA | Other |
| 29582 | (29582) 1998 FR58 | 42300 | (42300) 2001 UU140 | 2015-10-18 | 3.011e+06 | 0.044 | 8.922 | 3.895 | Other | MBA |
| 18210 | (18210) 4529 P-L | 57028 | (57028) 2000 UJ60 | 2014-07-25 | 4.165e+06 | 0.045 | 5.134 | 4.472 | Other | MBA |
| 4818 | (4818) Elgar | 80909 | (80909) 2000 DL59 | 2016-01-22 | 6.266e+06 | 0.046 | 6.768 | 1.864 | Other | Other |
| 39686 | (39686) Takeshihara | 153484 | (153484) 2001 RQ80 | 2014-07-24 | 3.280e+06 | 0.047 | 5.629 | 2.347 | MBA | Other |
| 24560 | (24560) 4517 P-L | 65553 | (65553) 1297 T-2 | 2016-01-18 | 5.774e+06 | 0.047 | 6.463 | 2.347 | Other | MBA |
| 11298 | (11298) Gide | 202763 | (202763) 2007 RW94 | 2014-08-14 | 1.479e+06 | 0.048 | 5.895 | 2.241 | MBA | Other |
| 7762 | (7762) 1990 SY2 | 21027 | (21027) 1989 SR5 | 2016-06-11 | 6.396e+06 | 0.048 | 9.342 | 4.903 | MBA | Other |
| 86044 | (86044) 1999 OD2 | 444114 | (444114) 2004 TT132 | 2016-06-08 | 6.605e+06 | 0.050 | 5.134 | 2.347 | Other | Other |
| 3987 | (3987) Wujek | 216385 | (216385) 2008 BX49 | 2015-11-14 | 5.964e+06 | 0.050 | 13.504 | 1.550 | Other | Other |
| 57480 | (57480) 2001 SO153 | 141023 | (141023) 2001 WV51 | 2014-08-09 | 2.740e+06 | 0.050 | 5.134 | 2.241 | Other | Other |
| 8096 | (8096) Emilezola | 119107 | (119107) 2001 OX56 | 2017-01-25 | 5.969e+06 | 0.051 | 6.172 | 1.952 | Other | MBA |
| 49345 | (49345) 1998 WH4 | 359737 | (359737) 2011 UR50 | 2014-07-24 | 6.963e+06 | 0.051 | 7.087 | 1.414 | Other | Other |
| 52099 | (52099) 2589 P-L | 144615 | (144615) 2004 FD62 | 2014-09-26 | 6.933e+06 | 0.053 | 5.895 | 2.573 | Other | Other |
| 1039 | (1039) Sonneberga | 204753 | (204753) 2006 JK24 | 2016-07-29 | 6.688e+06 | 0.053 | 17.802 | 1.780 | Other | Other |
| 3485 | (3485) Barucci | 135643 | (135643) 2002 JM109 | 2015-11-18 | 6.963e+06 | 0.054 | 9.342 | 1.232 | Other | Other |
| 16115 | (16115) 1999 XH25 | 303927 | (303927) 2005 UV392 | 2016-11-08 | 7.477e+06 | 0.056 | 8.520 | 1.780 | Other | Other |
| 4389 | (4389) Durbin | 424808 | (424808) 2008 UB107 | 2016-11-19 | 2.922e+06 | 0.056 | 12.896 | 1.481 | Other | Other |
| 3019 | (3019) Kulin | 85569 | (85569) 1998 BG18 | 2015-09-16 | 5.157e+06 | 0.057 | 14.807 | 2.821 | MBA | Other |
| 77714 | (77714) 2001 OY47 | 82945 | (82945) 2001 QN117 | 2016-06-27 | 5.553e+06 | 0.057 | 5.895 | 4.270 | Other | Other |
| 15498 | (15498) 1999 EQ4 | 31734 | (31734) 1999 JT71 | 2017-05-28 | 4.487e+06 | 0.059 | 11.761 | 6.768 | MBA | MBA |
| 11170 | (11170) 1998 FY34 | 20528 | (20528) Kyleyawn | 2015-08-28 | 6.961e+06 | 0.059 | 6.463 | 4.903 | MBA | MBA |


### 2c. Grande + lento + cercano (máximo interés físico)

Al menos un cuerpo D ≳ 50 km, v_rel ≤ 1 km/s, dist ≤ 1e+06 km.

| number_1 | designation_1 | number_2 | designation_2 | date_utc | dist_km | rel_vel_km_s | diameter_1_km | diameter_2_km | class_1 | class_2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 135 | (135) Hertha | 281202 | (281202) 2007 FE47 | 2016-09-18 | 4.556e+05 | 0.579 | 80.253 | 1.024 | MBA | Other |
| 371 | (371) Bohemia | 14711 | (14711) 2000 CG36 | 2016-05-14 | 5.920e+05 | 0.982 | 64.041 | 5.895 | MBA | Other |
| 83 | (83) Beatrix | 170676 | (170676) 2003 YE180 | 2014-09-18 | 6.098e+05 | 0.684 | 65.836 | 2.140 | MBA | Other |
| 678 | (678) Fredegundis | 200373 | (200373) 2000 QH89 | 2015-03-13 | 6.271e+05 | 0.964 | 55.778 | 2.044 | Other | MBA |
| 9 | (9) Metis | 345254 | (345254) 2005 UT477 | 2016-10-30 | 7.138e+05 | 0.979 | 197 | 1.232 | MBA | Other |
| 639 | (639) Latona | 10793 | (10793) Quito | 2014-12-18 | 7.273e+05 | 0.604 | 81.369 | 10.727 | Other | Other |
| 135 | (135) Hertha | 93002 | (93002) 2000 RN85 | 2017-04-26 | 7.668e+05 | 0.953 | 80.253 | 1.414 | MBA | Other |
| 91 | (91) Aegina | 388085 | (388085) 2005 UA113 | 2015-05-06 | 7.770e+05 | 0.908 | 60.598 | 0.934 | Other | Other |
| 431 | (431) Nephele | 368859 | (368859) 2006 JX39 | 2015-09-11 | 7.947e+05 | 0.428 | 64.041 | 1.952 | MBA | Other |
| 77 | (77) Frigga | 1808 | (1808) Bellerophon | 2016-06-28 | 8.087e+05 | 0.826 | 70.220 | 13.504 | Other | MBA |
| 377 | (377) Campania | 435772 | (435772) 2008 UF285 | 2016-01-18 | 8.772e+05 | 0.825 | 59.219 | 1.290 | Other | Other |
| 363 | (363) Padua | 62571 | (62571) 2000 SY274 | 2016-09-14 | 9.212e+05 | 0.894 | 56.035 | 3.719 | MBA | MBA |
| 222 | (222) Lucia | 426307 | (426307) 2012 TP165 | 2015-10-17 | 9.649e+05 | 0.760 | 53.023 | 1.864 | Other | Other |
| 16 | (16) Psyche | 139485 | (139485) 2001 PT14 | 2015-01-19 | 9.661e+05 | 0.869 | 235 | 2.573 | MBA | MBA |
| 7 | (7) Iris | 153539 | (153539) 2001 SD101 | 2015-08-17 | 9.968e+05 | 0.806 | 281 | 2.347 | MBA | Other |


## 3. Pares en la misma región dinámica (proxy de familia)

Pares cuyos elementos osculantes `(a, e, i)` son mutuamente cercanos (Δa/a ≤ 1 %, Δe ≤ 0.02, Δi ≤ 1°) y que además tuvieron un encuentro físico cercano. Se añaden las columnas de elementos de ambos cuerpos.

| number_1 | designation_1 | number_2 | designation_2 | date_utc | dist_km | rel_vel_km_s | diameter_1_km | diameter_2_km | class_1 | class_2 | a1 | a2 | e1 | e2 | i1 | i2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 200764 | (200764) 2001 XP3 | 131597 | (131597) 2001 XT2 | 2014-08-13 | 49511 | 1.794 | 2.954 | 2.821 | MBA | MBA | 2.297 | 2.289 | 0.189 | 0.203 | 25.146 | 24.460 |
| 89127 | (89127) 2001 UX3 | 90256 | (90256) 2003 BQ72 | 2015-01-30 | 1.210e+05 | 4.588 | 3.392 | 2.694 | MBA | MBA | 2.444 | 2.444 | 0.099 | 0.087 | 6.898 | 7.235 |
| 41528 | (41528) 2000 RN4 | 62757 | (62757) 2000 UU8 | 2015-06-08 | 1.418e+05 | 4.996 | 3.239 | 2.573 | MBA | MBA | 2.403 | 2.413 | 0.154 | 0.142 | 7.098 | 6.756 |
| 110727 | (110727) 2001 TR235 | 280694 | (280694) 2005 FX4 | 2016-10-20 | 1.513e+05 | 5.218 | 4.682 | 3.239 | MBA | MBA | 3.151 | 3.132 | 0.184 | 0.185 | 11.860 | 10.875 |
| 134682 | (134682) 1999 XM27 | 136228 | (136228) 2003 WN106 | 2014-08-23 | 1.517e+05 | 4.560 | 2.457 | 2.044 | MBA | MBA | 2.370 | 2.354 | 0.210 | 0.194 | 6.067 | 6.231 |
| 81211 | (81211) 2000 FX16 | 270155 | (270155) 2001 SN106 | 2015-05-13 | 1.533e+05 | 5.943 | 3.552 | 2.044 | MBA | MBA | 2.361 | 2.359 | 0.266 | 0.275 | 8.123 | 8.730 |
| 92384 | (92384) 2000 HS75 | 194052 | (194052) 2001 SY106 | 2015-07-19 | 1.640e+05 | 5.293 | 2.140 | 1.864 | MBA | MBA | 2.319 | 2.309 | 0.090 | 0.080 | 6.606 | 6.116 |
| 36099 | (36099) 1999 RE113 | 265722 | (265722) 2005 UX383 | 2015-10-15 | 1.789e+05 | 4.956 | 4.270 | 2.457 | MBA | MBA | 2.648 | 2.628 | 0.146 | 0.132 | 12.998 | 12.043 |
| 64618 | (64618) 2001 XQ28 | 121048 | (121048) 1999 CF42 | 2016-04-16 | 1.816e+05 | 3.932 | 2.694 | 2.694 | MBA | MBA | 2.278 | 2.273 | 0.179 | 0.167 | 6.619 | 6.264 |
| 33376 | (33376) Medi | 97800 | (97800) 2000 NO23 | 2016-12-20 | 1.883e+05 | 3.962 | 3.552 | 2.694 | MBA | MBA | 2.215 | 2.218 | 0.120 | 0.103 | 6.278 | 6.666 |
| 76929 | (76929) 2001 AX34 | 159931 | (159931) 2005 VY5 | 2015-01-17 | 1.888e+05 | 6.657 | 4.682 | 1.700 | MBA | MBA | 1.931 | 1.916 | 0.070 | 0.076 | 19.599 | 19.527 |
| 47734 | (47734) 2000 DX55 | 203726 | (203726) 2002 QW35 | 2014-11-02 | 1.973e+05 | 1.863 | 3.719 | 2.694 | MBA | MBA | 2.787 | 2.761 | 0.071 | 0.072 | 4.002 | 4.217 |
| 80597 | (80597) 2000 AD147 | 188253 | (188253) 2002 XM73 | 2015-08-30 | 2.121e+05 | 0.868 | 2.694 | 1.864 | MBA | MBA | 2.262 | 2.283 | 0.038 | 0.052 | 6.292 | 6.057 |
| 28304 | (28304) 1999 CC75 | 138282 | (138282) 2000 GA30 | 2014-12-12 | 2.127e+05 | 7.094 | 8.520 | 3.552 | MBA | MBA | 2.685 | 2.708 | 0.153 | 0.141 | 10.963 | 11.414 |
| 67506 | (67506) 2000 RL47 | 72110 | (72110) 2000 YR55 | 2016-05-17 | 2.132e+05 | 5.205 | 2.694 | 2.457 | MBA | MBA | 2.323 | 2.318 | 0.139 | 0.154 | 6.490 | 5.553 |
| 67529 | (67529) 2000 RQ90 | 252066 | (252066) 2000 ST139 | 2015-11-09 | 2.157e+05 | 2.369 | 2.457 | 1.481 | MBA | MBA | 2.237 | 2.239 | 0.107 | 0.101 | 6.099 | 6.197 |
| 16271 | (16271) Duanenichols | 226318 | (226318) 2003 DG18 | 2017-05-09 | 2.243e+05 | 2.096 | 3.719 | 1.700 | MBA | MBA | 2.426 | 2.446 | 0.141 | 0.144 | 4.527 | 4.010 |
| 25790 | (25790) 2000 CW57 | 90491 | (90491) 2004 DW22 | 2017-01-02 | 2.271e+05 | 5.636 | 8.137 | 3.239 | MBA | MBA | 2.780 | 2.795 | 0.161 | 0.176 | 9.547 | 8.735 |
| 146272 | (146272) 2001 EY5 | 79277 | (79277) 1995 SB25 | 2016-07-06 | 2.321e+05 | 3.632 | 2.573 | 2.241 | MBA | MBA | 2.308 | 2.307 | 0.153 | 0.143 | 7.150 | 8.098 |
| 57205 | (57205) 2001 QM55 | 118514 | (118514) 2000 DJ104 | 2015-05-30 | 2.322e+05 | 5.662 | 2.347 | 2.044 | MBA | MBA | 2.282 | 2.291 | 0.126 | 0.131 | 6.929 | 5.929 |
| 79951 | (79951) 1999 CY93 | 84605 | (84605) 2002 VL34 | 2017-05-10 | 2.374e+05 | 1.208 | 3.239 | 3.239 | MBA | MBA | 2.720 | 2.712 | 0.098 | 0.097 | 6.188 | 5.913 |
| 50084 | (50084) 2000 AZ89 | 84486 | (84486) 2002 TY275 | 2017-05-23 | 2.384e+05 | 3.264 | 3.392 | 2.347 | MBA | MBA | 2.306 | 2.299 | 0.062 | 0.061 | 7.064 | 6.205 |
| 28115 | (28115) 1998 SN50 | 67188 | (67188) 2000 CV28 | 2016-08-18 | 2.401e+05 | 3.803 | 3.239 | 2.241 | MBA | MBA | 2.373 | 2.384 | 0.153 | 0.162 | 1.427 | 1.430 |
| 27320 | (27320) Vellinga | 186783 | (186783) 2004 DT48 | 2015-05-13 | 2.485e+05 | 2.389 | 2.694 | 1.624 | MBA | MBA | 2.388 | 2.371 | 0.133 | 0.129 | 0.635 | 0.864 |
| 145595 | (145595) 2006 QE3 | 115572 | (115572) 2003 UY85 | 2016-07-27 | 2.493e+05 | 3.079 | 2.241 | 2.140 | MBA | MBA | 2.216 | 2.211 | 0.186 | 0.192 | 8.157 | 8.294 |
| 113914 | (113914) 2002 TF284 | 142558 | (142558) 2002 TP61 | 2014-08-14 | 2.512e+05 | 0.913 | 2.241 | 1.780 | MBA | MBA | 2.307 | 2.285 | 0.107 | 0.116 | 6.285 | 6.423 |
| 52013 | (52013) 2002 LJ59 | 244498 | (244498) 2002 TF84 | 2014-09-07 | 2.515e+05 | 7.681 | 4.903 | 2.694 | MBA | MBA | 2.766 | 2.790 | 0.205 | 0.203 | 15.298 | 16.050 |
| 68100 | (68100) 2000 YV120 | 141631 | (141631) 2002 JL49 | 2015-01-11 | 2.542e+05 | 3.478 | 2.241 | 1.952 | MBA | MBA | 2.271 | 2.281 | 0.151 | 0.138 | 4.416 | 3.948 |
| 29453 | (29453) 1997 RU6 | 125513 | (125513) 2001 WN39 | 2016-02-17 | 2.575e+05 | 5.756 | 3.719 | 2.347 | MBA | MBA | 2.401 | 2.416 | 0.149 | 0.161 | 5.558 | 5.770 |
| 27320 | (27320) Vellinga | 79650 | (79650) 1998 SB16 | 2015-03-31 | 2.655e+05 | 3.222 | 2.694 | 2.044 | MBA | MBA | 2.388 | 2.390 | 0.133 | 0.130 | 0.635 | 0.429 |


## 4. Perturbadores

### 4a. Encuentros que tocan uno de los 16 perturbadores estudiados

Lista de referencia (ya cubierta): 1, 2, 3, 4, 7, 10, 15, 16, 31, 52, 65, 87, 88, 107, 511, 704. Sirve para **separar lo ya trabajado** del descubrimiento.

| number_1 | designation_1 | number_2 | designation_2 | date_utc | dist_km | rel_vel_km_s | diameter_1_km | diameter_2_km | class_1 | class_2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 65 | (65) Cybele | 250012 | (250012) 2002 AX51 | 2016-05-20 | 1.334e+05 | 2.525 | 168 | 2.347 | MBA | Other |
| 4 | (4) Vesta | 446361 | (446361) 2014 HH24 | 2014-11-03 | 1.401e+05 | 4.159 | 814 | 1.073 | MBA | Other |
| 4 | (4) Vesta | 387639 | (387639) 2002 QL8 | 2016-10-27 | 2.013e+05 | 1.999 | 814 | 1.073 | MBA | Other |
| 52 | (52) Europa | 88978 | (88978) 2001 TR62 | 2015-09-14 | 2.846e+05 | 2.988 | 194 | 3.392 | Other | Other |
| 4 | (4) Vesta | 310169 | (310169) 2011 SD7 | 2015-01-04 | 3.145e+05 | 5.862 | 814 | 1.952 | MBA | Other |
| 31 | (31) Euphrosyne | 382982 | (382982) 2005 CR15 | 2017-01-28 | 3.656e+05 | 9.043 | 159 | 1.550 | Other | Other |
| 15 | (15) Eunomia | 5199 | (5199) Dortmund | 2015-06-16 | 4.080e+05 | 5.688 | 312 | 14.140 | Other | MBA |
| 88 | (88) Thisbe | 400715 | (400715) 2009 SM15 | 2014-11-04 | 4.117e+05 | 2.439 | 139 | 1.624 | MBA | Other |
| 704 | (704) Interamnia | 425470 | (425470) 2010 EG107 | 2016-09-03 | 4.419e+05 | 7.686 | 230 | 1.481 | MBA | Other |
| 3 | (3) Juno | 278396 | (278396) 2007 PD47 | 2016-01-03 | 4.682e+05 | 4.658 | 305 | 2.044 | Other | Other |
| 704 | (704) Interamnia | 432371 | (432371) 2009 WH87 | 2016-07-22 | 4.720e+05 | 5.565 | 230 | 1.073 | MBA | Other |
| 87 | (87) Sylvia | 55420 | (55420) 2001 TV20 | 2016-04-09 | 4.750e+05 | 8.885 | 145 | 8.922 | MBA | MBA |
| 15 | (15) Eunomia | 284181 | (284181) 2006 AT56 | 2016-09-11 | 4.775e+05 | 5.034 | 312 | 1.073 | Other | Other |
| 4 | (4) Vesta | 450461 | (450461) 2005 WS15 | 2014-11-09 | 4.927e+05 | 5.798 | 814 | 0.777 | MBA | Other |
| 88 | (88) Thisbe | 454516 | (454516) 2014 ON214 | 2016-06-13 | 5.169e+05 | 7.377 | 139 | 2.140 | MBA | Other |
| 15 | (15) Eunomia | 24066 | (24066) Eriksorensen | 2014-12-12 | 5.288e+05 | 5.095 | 312 | 3.552 | Other | Other |
| 52 | (52) Europa | 273122 | (273122) 2006 FL54 | 2016-08-12 | 5.688e+05 | 5.103 | 194 | 2.241 | Other | Other |
| 1 | (1) Ceres | 35385 | (35385) 1997 WL37 | 2014-10-28 | 5.714e+05 | 4.600 | 763 | 2.044 | MBA | Other |
| 4 | (4) Vesta | 200427 | (200427) 2000 SA336 | 2017-02-28 | 5.943e+05 | 5.618 | 814 | 1.780 | MBA | MBA |
| 1 | (1) Ceres | 24836 | (24836) 1995 TO1 | 2015-03-01 | 5.977e+05 | 3.614 | 763 | 2.954 | MBA | Other |
| 15 | (15) Eunomia | 307407 | (307407) 2002 TP159 | 2016-07-20 | 5.989e+05 | 8.044 | 312 | 2.140 | Other | Other |
| 1 | (1) Ceres | 147856 | (147856) 2005 UH343 | 2017-01-08 | 6.018e+05 | 4.635 | 763 | 4.903 | MBA | MBA |
| 88 | (88) Thisbe | 26094 | (26094) 1988 NU | 2016-10-13 | 6.209e+05 | 1.282 | 139 | 5.629 | MBA | MBA |
| 4 | (4) Vesta | 125989 | (125989) 2001 YA31 | 2015-03-27 | 6.753e+05 | 2.043 | 814 | 1.780 | MBA | MBA |
| 1 | (1) Ceres | 65095 | (65095) 2002 CN3 | 2016-01-10 | 6.795e+05 | 5.567 | 763 | 5.895 | MBA | Other |
| 15 | (15) Eunomia | 452522 | (452522) 2004 TO40 | 2014-10-24 | 6.823e+05 | 3.874 | 312 | 0.934 | Other | Other |
| 16 | (16) Psyche | 288581 | (288581) 2004 HV65 | 2014-10-16 | 7.029e+05 | 3.092 | 235 | 1.290 | MBA | Other |
| 15 | (15) Eunomia | 310561 | (310561) 2001 QQ230 | 2016-07-30 | 7.048e+05 | 1.772 | 312 | 1.481 | Other | Other |
| 31 | (31) Euphrosyne | 42591 | (42591) 1997 GE42 | 2017-03-21 | 7.093e+05 | 10.114 | 159 | 3.552 | Other | MBA |
| 704 | (704) Interamnia | 236960 | (236960) 2007 UX74 | 2016-12-12 | 7.123e+05 | 5.916 | 230 | 2.140 | MBA | Other |


### 4b. Cuerpos grandes FUERA de los 16 (candidatos a masa nueva → F4)

Cuerpos con D ≳ 100 km, no incluidos en los 16, rankeados por número de encuentros *útiles* (v_rel ≤ 3 km/s y dist ≤ 3×10⁶ km). Más eventos buenos ⇒ mejor candidato a determinación de masa.

| big_number | big_name | big_diam_km | big_class | n_useful_events | best_dist_km | min_vrel_km_s |
| --- | --- | --- | --- | --- | --- | --- |
| 9 | (9) Metis | 197 | MBA | 37 | 7.138e+05 | 0.378 |
| 30 | (30) Urania | 109 | Other | 36 | 4.398e+05 | 0.817 |
| 40 | (40) Harmonia | 141 | Other | 36 | 6.578e+05 | 1.028 |
| 19 | (19) Fortuna | 133 | MBA | 32 | 7.433e+05 | 0.678 |
| 21 | (21) Lutetia | 120 | Other | 30 | 8.179e+05 | 0.609 |
| 64 | (64) Angelina | 104 | MBA | 30 | 8.731e+05 | 0.992 |
| 44 | (44) Nysa | 140 | Other | 30 | 8.946e+05 | 0.961 |
| 20 | (20) Massalia | 178 | Other | 28 | 8.778e+05 | 0.514 |
| 29 | (29) Amphitrite | 240 | Other | 27 | 3.166e+05 | 0.643 |
| 27 | (27) Euterpe | 141 | Other | 26 | 4.634e+05 | 0.847 |
| 26 | (26) Proserpina | 112 | MBA | 25 | 68766 | 0.834 |
| 128 | (128) Nemesis | 113 | MBA | 25 | 4.138e+05 | 1.036 |
| 11 | (11) Parthenope | 174 | MBA | 24 | 6.902e+05 | 0.627 |
| 103 | (103) Hera | 104 | Other | 23 | 8.434e+05 | 1.083 |
| 32 | (32) Pomona | 109 | Other | 21 | 4.301e+05 | 1.359 |
| 8 | (8) Flora | 179 | Other | 21 | 8.824e+05 | 0.525 |
| 63 | (63) Ausonia | 110 | Other | 18 | 9.372e+05 | 0.784 |
| 55 | (55) Pandora | 102 | Other | 17 | 3.203e+05 | 1.532 |
| 45 | (45) Eugenia | 114 | Other | 16 | 9.817e+05 | 1.016 |
| 106 | (106) Dione | 117 | Other | 15 | 4.441e+05 | 1.207 |
| 346 | (346) Hermentaria | 133 | Other | 15 | 4.958e+05 | 1.634 |
| 349 | (349) Dembowska | 232 | Other | 14 | 5.211e+05 | 1.172 |
| 192 | (192) Nausikaa | 133 | MBA | 14 | 7.727e+05 | 1.786 |
| 221 | (221) Eos | 104 | Other | 14 | 9.024e+05 | 1.047 |
| 5 | (5) Astraea | 152 | MBA | 14 | 1.046e+06 | 1.576 |
| 120 | (120) Lachesis | 100 | Other | 13 | 1.861e+05 | 0.845 |
| 68 | (68) Leto | 156 | MBA | 13 | 6.305e+05 | 1.213 |
| 37 | (37) Fides | 124 | Other | 13 | 6.395e+05 | 1.890 |
| 241 | (241) Germania | 108 | MBA | 12 | 1.065e+06 | 1.593 |
| 42 | (42) Isis | 111 | MBA | 11 | 6.760e+05 | 1.020 |


## 5. Candidatos a seguimiento

### 5a. Eventos genuinamente notables / potencialmente no catalogados

- **2 encuentros grande-grande D ≳ 100 km** en el catálogo (sección 1.2). Cualquier par de este grupo que **no** toque a los 16 perturbadores conocidos es un evento de alto perfil sin cobertura previa obvia; ver la marca de perturbador en §4a.
- El encuentro grande-grande (D ≳ 50 km) más cercano **fuera de los 16** es (305) Gordonia × (830) Petropolitana el 2014-09-18 a 9.876e+05 km — candidato a revisión N-cuerpos.
- Los encuentros más cercanos en términos absolutos (§2a) merecen revalidación N-cuerpos: a esas distancias la aproximación Kepler es más frágil, pero si sobreviven son los eventos geométricamente más notables del dataset.

### 5b. Cuerpos grandes fuera de los 16 con más encuentros útiles (→ F4)

Ranking del §4b (top 10). Estos son los perturbadores no estudiados con más eventos de baja v_rel y corta distancia, es decir, con el mayor potencial de deflexión medible para una masa nueva:

- **(9) Metis** (D≈197 km, clase MBA): 37 eventos útiles; mejor par a 7.138e+05 km, v_rel mín 0.378 km/s.
- **(30) Urania** (D≈109 km, clase Other): 36 eventos útiles; mejor par a 4.398e+05 km, v_rel mín 0.817 km/s.
- **(40) Harmonia** (D≈141 km, clase Other): 36 eventos útiles; mejor par a 6.578e+05 km, v_rel mín 1.028 km/s.
- **(19) Fortuna** (D≈133 km, clase MBA): 32 eventos útiles; mejor par a 7.433e+05 km, v_rel mín 0.678 km/s.
- **(21) Lutetia** (D≈120 km, clase Other): 30 eventos útiles; mejor par a 8.179e+05 km, v_rel mín 0.609 km/s.
- **(64) Angelina** (D≈104 km, clase MBA): 30 eventos útiles; mejor par a 8.731e+05 km, v_rel mín 0.992 km/s.
- **(44) Nysa** (D≈140 km, clase Other): 30 eventos útiles; mejor par a 8.946e+05 km, v_rel mín 0.961 km/s.
- **(20) Massalia** (D≈178 km, clase Other): 28 eventos útiles; mejor par a 8.778e+05 km, v_rel mín 0.514 km/s.
- **(29) Amphitrite** (D≈240 km, clase Other): 27 eventos útiles; mejor par a 3.166e+05 km, v_rel mín 0.643 km/s.
- **(27) Euterpe** (D≈141 km, clase Other): 26 eventos útiles; mejor par a 4.634e+05 km, v_rel mín 0.847 km/s.

> Contrastar contra los candidatos F4 propuestos en `docs/mass_layer_f4_design.md` (24 Themis, 532 Herculina, 29 Amphitrite, 354 Eleonora) y priorizar los que además aparezcan alto en este ranking.

---

_Generado por `scripts/bench/mine_notable_encounters.py` el 2026-07-02._
