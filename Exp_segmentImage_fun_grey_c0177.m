function [BW,maskedImage] = segmentImage(X)
%segmentImage Segment image using auto-generated code from imageSegmenter app
%  [BW,MASKEDIMAGE] = segmentImage(X) segments image X using auto-generated
%  code from the imageSegmenter app. The final segmentation is returned in
%  BW, and a masked image is returned in MASKEDIMAGE.

% Auto-generated by imageSegmenter app on 24-Jul-2022
%----------------------------------------------------


% Adjust data to span data range.
X = imadjust(X);

% Graph cut
foregroundInd = [563611 563614 563615 563616 563619 564685 564687 564688 564701 564702 566842 566844 566864 567921 567945 570078 570080 570107 571158 571190 572238 574397 575477 575511 576595 577635 580875 580918 583034 583080 586274 586324 588433 588485 589513 591727 592753 593831 596051 597132 598151 600374 601390 602470 602534 605708 605775 606857 606858 607868 607938 610099 611106 613343 614346 616505 618744 619826 620825 623067 624063 626309 627303 628470 629462 631712 632702 632793 633782 637020 637115 639180 640360 641373 641376 641377 641380 641382 642420 642448 642450 642452 642520 643499 643545 643548 643550 643601 644605 644606 644633 644634 644682 645683 646739 646797 646842 647879 647924 648897 650002 650043 650085 651126 652137 653216 653240 653290 653327 655399 655451 655488 657536 657613 658696 658730 659717 660773 660891 661974 662933 662955 663017 664035 666172 666261 666296 667274 667376 667377 669412 669434 669539 670514 670584 670619 670620 671700 673754 673826 673860 673862 674810 674942 674943 675912 678050 678148 679152 683450 685632 686690 686791 687770 689952 694250 694271 695351 697488 698568 698591 699757 701918 705048 705071 708311 709391 709480 712608 713688 714882 716928 718029 720168 725568 725589 726648 726669 730004 730968 736368 736387 739608 740707 740804 745027 746088 747285 748247 748267 751487 753665 754845 756885 756905 759065 761204 762305 766625 767807 773104 773207 777424 779689 780662 782929 784009 784982 789411 790382 792651 794700 794811 798052 799019 802372 803339 804417 805612 807772 811976 815216 818453 818574 819654 820613 822894 827091 828171 830330 832616 833570 837888 840048 841256 843287 849765 850977 853003 856377 858537 861640 862719 865019 867179 870277 871356 872581 876754 876901 882152 882302 883231 883384 886624 889865 891869 892027 895267 897265 898509 902831 903743 903911 910219 912555 915795 917956 924439 925333 930921 933969 934161 935243 939365 942803 943683 947124 951240 953604 957716 959875 962246 966566 968511 971966 978226 980606 981463 988166 989021 991406 992258 993336 996806 1001126 1001972 1006291 1006526 1008447 1009766 1015166 1018162 1020566 1025717 1025966 1027875 1031112 1033526 1037589 1040826 1041904 1049461 1051619 1054042 1058095 1060522 1062414 1068888 1074284 1077800 1078603 1083200 1083996 1093712 1100478 1104504 1108821 1109118 1114217 1114518 1115598 1118534 1120998 1121772 1126398 1127169 1130718 1131484 1131797 1135037 1135803 1136879 1140437 1141197 1143677 1144436 1146592 1146917 1149077 1149831 1151987 1155225 1155557 1159542 1166017 1166357 1169597 1170333 1173571 1174995 1175727 1178965 1183282 1186520 1188677 1189035 1191915 1196232 1197673 1201627 1203786 1207023 1208101 1210631 1212417 1215655 1216031 1219971 1221048 1225366 1226828 1227523 1230759 1235078 1240471 1242629 1243026 1245866 1246944 1257735 1258811 1260304 1260970 1269596 1272182 1272833 1274992 1283621 1284060 1284698 1286855 1292248 1294404 1295481 1299176 1299795 1304576 1305185 1309500 1312736 1314294 1318131 1319210 1319692 1320771 1325678 1326169 1328913 1329991 1330487 1331567 1332144 1333223 1338615 1339124 1340767 1341845 1344002 1347759 1349396 1352630 1354234 1354785 1355863 1358017 1360174 1362870 1363406 1363948 1365562 1368797 1370951 1373108 1375261 1377981 1378492 1380649 1383884 1385537 1386040 1387115 1390934 1391425 1393579 1395252 1395733 1396811 1398964 1400650 1401118 1403274 1404348 1406502 1406505 1407126 1408658 1410365 1411889 1411890 1414044 1414045 1415116 1416840 1417275 1420512 1421585 1423318 1423743 1424820 1427636 1429131 1429134 1431288 1432364 1433033 1434522 1435598 1435599 1437351 1437754 1442066 1442069 1442749 1443144 1447457 1448536 1450691 1451383 1451768 1457162 1457860 1459319 1462556 1464336 1464714 1467575 1467951 1469731 1470107 1472968 1473346 1475502 1479445 1482683 1483056 1484842 1488078 1488451 1495633 1496005 1499243 1503187 1503560 1509659 1510034 1510736 1515053 1516511 1518669 1521905 1522601 1526223 1526916 1530541 1532310 1532698 1535938 1537017 1538776 1541335 1542010 1545653 1546322 1548477 1550631 1551050 1552788 1554941 1557528 1560328 1562482 1562928 1564635 1567868 1570486 1572183 1574335 1575885 1577568 1579721 1580795 1582949 1584525 1589416 1589923 1591568 1595882 1598036 1600721 1604507 1606662 1608281 1608818 1612053 1614208 1618523 1619079 1619599 1621758 1624477 1624991 1629305 1633117 1635773 1639597 1641167 1642242 1646556 1648237 1651477 1655183 1657338 1659037 1661652 1662277 1664437 1664887 1665517 1667039 1668758 1670920 1671353 1673504 1673507 1677403 1677811 1679965 1680645 1681726 1682124 1684968 1685347 1685350 1687502 1690372 1691813 1692879 1692887 1696854 1697188 1699340 1700098 1701495 1703646 1704420 1705801 1706583 1706875 1709022 1709027 1709827 1711989 1712251 1714150 1714405 1715233 1716555 1717395 1718478 1718704 1720641 1720850 1720856 1722805 1722999 1723005 1723888 1725141 1725144 1725150 1725152 1727279 1727284 1727293 1727296 1728214 1728340 1728346 1728354 1730492 1730497 1732539 1732633 1732639 1732643 1734706 1734785 1735789 1735855 1735857 1735861 1738006 1738010 1739082 1739085 1741240 1742317 1744446 1744473 1744474 1744476 1745552 1745553 1746611 1748790 1748792 1749870 1752018 1753109 1754181 1755267 1756346 1757426 1758509 1761754 ];
backgroundInd = [32 33 35 38 41 42 45 49 53 56 57 60 64 65 66 69 73 76 77 80 207 208 210 768 907 1109 1110 1501 1531 2244 2246 2249 2252 2254 2257 2260 2263 2268 2271 2274 2371 2579 2612 2927 3268 3360 3363 3364 3367 3370 3372 3375 3440 3445 3454 3456 3925 3928 3931 3936 3944 3999 4003 4344 4778 5538 5541 5543 5589 5590 5597 5619 6077 6081 6105 6110 6113 6116 6120 6124 6132 6141 6142 6147 6150 6154 6158 6311 6502 6503 6629 6631 6635 6636 6639 6642 6645 6647 6650 6653 6654 6661 6701 6898 6945 7580 7974 8026 8229 8230 8233 8472 8864 9301 9304 9306 9554 9736 9739 9946 10189 10815 10816 11456 11460 11895 12109 12287 12531 12532 12795 13438 14053 14273 14276 14684 14688 15133 15357 15601 16594 16682 16841 17520 17522 17915 17917 17920 18373 18603 18606 18607 18610 18612 18623 18625 18628 18634 18640 18643 18651 18660 18663 18667 18673 18676 18682 18688 18703 18710 18722 18732 18742 19283 19314 19318 19319 19325 19326 19329 19331 19334 19927 20057 20062 20063 20071 20072 20533 21112 21117 21120 21129 21134 21471 21613 22088 22092 23174 23177 23178 23254 23260 23263 23264 23607 23608 23625 23628 23657 23771 23772 24703 24849 25341 25345 25348 25349 25352 25355 25359 25362 25363 25370 25372 25375 25378 25381 25383 25389 25392 25394 25403 25406 25408 25772 25775 25821 25928 26857 26860 26904 28083 28085 29066 29067 29162 30150 30152 30155 31322 32317 32320 32321 32323 32324 32402 33407 33411 35642 36654 36722 38817 38882 39059 39062 39065 39069 39072 39079 39962 40983 41042 41247 41248 41251 41257 41262 41265 41267 41270 41274 41280 41282 41505 41513 41529 41535 41539 41544 41547 41549 41552 41557 41563 41566 42366 42368 42371 42375 43202 43659 43732 43740 44229 44538 44542 44545 44546 44548 44551 44554 44560 44562 44566 44572 44575 44578 44583 44586 44587 44590 44591 44594 45808 45911 45920 46398 46400 46403 46404 46405 46442 46758 46762 46765 46768 46770 46773 46884 47522 47854 47861 48086 48087 48565 50024 50025 50032 50038 51124 51125 51128 51130 51192 51195 51342 51805 52267 52922 53296 53301 53507 54003 54385 54388 54392 54398 54401 54406 54414 54417 55045 55678 57243 57844 57852 57858 57861 58285 59367 59405 60028 60035 60039 60079 60084 60087 60091 60095 60099 60102 60105 60110 60114 60115 60125 60131 60133 60136 60139 60144 60147 60148 60154 60244 60245 60248 60249 60251 60254 60257 60260 60269 60271 60272 60274 61527 62204 62211 62215 62218 62224 62227 62231 62235 62236 62318 62319 62327 62331 62339 62341 62344 62355 62365 62376 62385 62389 62396 62645 63725 64768 68008 68045 72365 75606 78846 80966 81006 84248 85328 88568 89605 89751 89969 89976 89987 89997 90004 90006 90011 90018 90020 90028 90032 91114 91128 91133 91808 91905 92092 92097 92115 92123 93307 93317 93332 93968 94104 94119 94136 94157 94188 94199 94216 94224 94763 94781 94786 94802 94819 94828 94838 94859 94870 95142 95504 95509 95519 96083 96223 96236 96247 96904 96907 96912 97052 97210 97691 97697 97712 97717 97728 97752 97767 97787 97799 97827 97834 97844 97866 97871 97891 97898 97903 97908 98378 99045 99056 99459 100087 100093 100098 100103 100109 100116 100117 100450 100537 101483 101617 102469 102610 103776 104856 105851 106804 106809 108011 108096 108970 111137 111251 111336 112218 112281 112416 113411 115463 115521 115656 116651 117628 118711 119891 120873 120971 121056 123081 124117 124211 126278 126321 126371 126454 129559 129611 132851 132934 134924 135093 136091 138164 138251 138333 139279 139331 140359 140411 141573 142571 142653 144731 146839 146891 147971 148053 148965 150079 150131 151159 151293 152239 153371 154401 154451 154533 157607 157644 157691 157773 158771 161013 162011 163173 164126 164171 164253 165206 166413 167409 168447 169489 170649 170733 171813 172769 172809 176011 176049 177049 177092 180369 180453 181412 181449 181533 183574 184773 185769 189009 190015 190057 190089 190173 190721 190724 190729 190743 190751 190761 190769 192217 192875 192941 192949 192957 192962 192968 192979 192988 192991 192996 193005 193015 193025 193035 193043 193045 193051 193053 193056 193061 193329 194149 194153 195019 195457 195573 196080 196086 196315 196321 196497 196569 197038 197054 197062 197067 197083 197089 197403 198471 198474 198480 198483 198697 198729 199180 199185 199255 199258 199263 199268 199304 200349 200352 200355 200361 200371 200621 200973 201937 201969 203470 203478 203483 204923 204928 205209 205616 205625 206001 208051 208056 208061 208069 208078 208089 208093 208097 208152 208417 208449 208845 210204 210261 210268 210279 210287 210293 210300 210538 210981 210993 210997 211657 212358 212768 212851 213121 213126 213134 213436 215592 215938 217088 218168 218511 219331 219585 220295 221408 222058 222061 223568 226775 226808 227888 228212 228524 228898 230131 231128 232138 232838 234368 236067 236528 237538 238223 239698 239768 240377 240777 240813 240848 242234 242529 244015 244088 244171 244684 245758 246213 246248 247255 249488 250568 251731 252224 254097 254888 255459 255968 257613 258052 258688 259208 259291 261331 261368 261645 262531 264530 265153 265960 266731 268535 268540 268555 268891 269464 270008 270686 270734 271617 272251 272899 273166 273248 273508 274290 275000 275077 275084 275491 275925 276488 276571 277249 277255 277530 278340 278564 278610 279312 279979 280770 281850 281889 282388 283746 284010 284131 284290 285129 286447 287200 287934 288330 288369 288451 289169 291011 291570 292689 293324 294001 294753 294811 296011 297234 297480 298089 299251 299999 300214 300872 301329 302629 303379 303454 304569 305417 305868 306696 306729 307177 307334 307891 308769 309971 310051 310846 311018 311051 311649 313028 313043 313054 313062 313069 313422 315341 315800 315963 316451 319662 319898 320851 322932 324012 324589 328305 328332 328741 328903 329614 330494 333056 335971 336945 336974 337368 338054 340185 340763 341294 341371 341491 342761 346665 346694 348311 349934 352065 352625 353174 353251 354627 355334 356607 358943 359097 360705 361785 362894 364338 366134 366651 368294 368371 368808 370425 370809 371534 372614 373771 375823 375854 378362 378517 379094 381254 384494 384993 385543 387734 387811 388072 390388 390943 392054 393211 393625 395482 396343 396374 399943 401774 402823 403181 403934 404416 404748 404755 404756 404759 404762 404764 404767 404771 404774 404779 404783 405864 405870 406063 407143 407174 407980 407983 408033 408331 408732 409058 409118 410383 410414 411215 411281 411284 412292 412367 412543 413654 414128 414529 415527 415783 415814 416692 417684 417777 419023 419052 419369 419840 419940 420447 420919 421021 422263 422292 422371 424153 425532 425842 426310 426428 427082 427389 427663 427692 427771 429546 430158 430624 430752 430903 430932 431703 432478 432915 433861 434143 434172 434938 435223 436332 437237 437385 437599 438177 438318 439255 440625 440652 440950 441414 441560 442491 442493 442891 443721 443892 444650 444947 444972 445728 445730 446805 446807 448044 448210 448291 448673 448963 448964 450663 451427 451450 452366 454069 454773 455607 455748 455959 456504 456690 456850 457930 458658 458988 460068 460376 462091 462250 464692 465490 467495 467630 468710 470088 472053 472896 473049 474112 475209 475399 475660 476138 477641 478301 478449 481543 481689 481958 482137 482769 484783 485865 486009 486094 487072 490186 490414 490522 491409 492489 494508 494903 495589 495712 495729 495816 498553 498969 499911 500712 501113 501378 503154 503378 503950 504369 504458 505449 505887 507186 507476 507593 511798 511913 512583 513443 513961 514180 515169 515402 515821 517204 517313 517980 518409 520795 521525 521846 523689 523809 524239 524455 526060 526615 526930 527268 528012 528129 530490 530932 531255 531353 532661 533416 533622 533733 534171 534815 535689 535894 535895 536330 536658 536753 536769 536970 536971 536973 536979 537408 538822 538913 539127 539128 539567 540647 541725 542063 542169 543884 546473 546582 548202 548633 549630 549729 550145 551889 552519 552874 554033 555038 555757 556122 556620 557273 557914 558353 558369 560443 561152 561527 561592 563688 563752 564390 564771 564827 564831 564851 566932 566981 568053 569096 569784 570179 570183 570185 570189 570196 570199 570204 570210 571331 572411 572504 573892 574103 575181 578418 578891 582131 585371 586846 587053 589209 592448 594011 596764 597251 600492 601079 601666 602652 604117 605398 607554 609132 610790 612372 612753 614532 615108 616185 619932 620310 623172 623738 626974 628666 629133 630732 631812 632182 633448 635605 638292 638843 641000 641532 643692 643786 644053 644238 646394 646932 648012 651252 652691 652870 657732 660972 661503 662052 662581 663226 663484 665820 666372 668878 670136 670692 673197 673932 674026 676610 679332 679849 682572 685066 685243 687972 688481 689383 693372 694956 699852 700173 703092 703590 704266 706648 708492 714383 715280 716052 717439 717621 721452 724095 724692 724993 730092 730569 731266 732547 734412 735784 737652 738732 741361 741972 742066 746292 746757 748732 749532 750612 751074 753852 757092 757548 758445 759252 763572 764919 766812 767262 769236 770052 771132 774817 775454 776534 776626 778054 779133 779774 781106 784530 785174 785422 787768 788414 789494 789925 793163 793816 797056 797295 797480 798558 801376 802456 802692 802875 805033 805696 806864 808089 808936 811098 812587 814338 815418 816901 817803 819058 821898 823374 824058 826439 826612 827298 828378 828769 832005 834164 837018 837234 838481 840258 841338 841421 842798 843710 847114 847818 849106 849978 850351 853218 854298 854666 856541 857742 862220 863298 864018 866535 867258 868338 875898 877171 877328 879486 882375 882725 884729 884881 886779 888120 888855 892095 892436 893512 896415 897830 899655 901066 903079 907215 909375 909557 909701 913777 916175 916935 918015 918193 921255 922335 922652 925657 925890 926655 932055 934524 935294 936457 937621 939614 942854 943018 943934 948254 948553 951494 953737 954734 955970 956894 957189 960134 961214 962377 963663 965534 965684 966762 969058 972014 972157 973177 975254 975534 978494 980654 982814 983088 983894 984029 986054 986325 987134 989294 992617 993878 997934 998057 998197 999275 1000214 1001174 1001434 1004671 1005582 1006688 1007652 1007765 1007909 1008824 1009923 1010892 1011000 1012225 1013146 1013158 1014235 1016389 1016390 1016392 1016393 1017620 1021692 1023017 1024932 1026255 1028172 1031649 1033572 1034886 1036812 1037044 1038972 1042212 1043519 1045452 1045678 1048916 1052154 1053012 1055393 1057332 1058412 1058631 1061652 1065972 1066186 1069212 1071582 1072452 1073532 1076772 1078932 1080012 1081296 1085614 1086490 1088650 1092091 1092970 1094050 1097290 1099450 1100727 1102690 1105930 1107010 1108283 1111522 1115650 1116918 1118890 1122130 1124477 1125370 1129874 1130770 1134010 1135090 1136350 1138330 1140490 1141747 1143730 1144984 1146970 1148050 1150381 1151290 1153450 1154699 1155777 1156690 1159016 1159930 1163170 1165493 1166571 1167490 1169650 1170730 1170890 1171968 1173970 1175207 1177210 1179370 1179523 1182610 1183690 1183839 1185850 1188155 1189090 1189234 1190170 1192470 1193410 1196650 1196788 1198810 1203130 1203261 1206370 1207450 1207579 1209610 1210690 1210816 1213930 1217170 1219330 1219452 1220410 1222570 1223650 1224847 1225810 1226890 1227970 1228085 1230244 1231210 1233370 1233482 1234561 1237690 1238877 1240928 1243088 1244274 1246328 1247408 1248488 1249671 1250648 1252909 1253888 1253988 1254968 1257226 1259288 1259385 1262528 1262623 1265768 1270178 1271167 1271256 1272247 1275487 1275575 1276567 1277647 1279807 1279893 1283047 1283129 1284127 1285288 1288447 1291687 1291764 1293922 1299247 1300399 1302558 1303636 1304647 1306876 1307887 1309035 1312207 1313287 1313353 1316590 1317607 1320847 1320908 1321987 1325167 1327327 1327384 1330567 1330624 1331702 1333807 1334887 1337101 1340285 1340339 1343525 1344659 1346817 1347845 1348925 1353245 1354374 1355454 1357565 1358693 1360850 1364045 1364090 1367285 1367330 1368409 1369445 1371649 1372685 1373765 1375967 1377005 1379165 1380286 1383485 1384606 1388885 1388926 1389967 1390004 1393244 1396484 1398607 1398644 1401847 1405122 1406202 1409406 1410489 1411566 1412646 1413729 1413762 1414841 1415886 1418046 1420206 1420209 1421286 1421321 1424561 1425606 1425640 1426686 1428880 1429926 1431006 1432120 1435326 1436438 1438566 1439646 1439678 1440726 1443996 1445046 1446156 1448286 1449366 1451526 1452606 1452635 1453686 1454795 1458006 1459113 1460166 1460193 1463432 1465592 1468805 1468830 1472045 1472068 1473125 1475308 1477443 1477467 1480683 1480705 1481784 1483923 1483944 1485024 1487163 1487183 1488243 1488261 1490403 1492563 1492580 1494723 1494738 1495804 1495818 1496884 1496896 1499044 1499054 1500124 1500131 1500133 1502284 1502286 1502287 1502289 1502290 1503369 1506609 1510929 1517407 1519567 1521692 1521695 1521696 1521698 1521699 1522807 1523846 1523859 1523887 1524939 1527099 1527127 1530321 1530367 1531419 1532478 1532499 1533607 1534687 1535739 1536794 1537899 1537927 1541167 1542219 1544407 1546506 1546541 1546567 1547647 1549806 1550861 1551901 1554101 1554126 1556219 1558422 1559455 1559526 1560582 1562693 1562766 1563825 1566006 1567011 1569170 1569246 1570305 1571385 1572408 1572486 1573487 1573566 1576804 1578885 1578964 1580026 1582204 1583201 1585444 1586440 1588666 1590842 1591922 1593002 1596152 1596242 1597230 1598388 1600562 1604868 1604881 1605867 1605961 1608108 1609201 1611261 1611361 1612441 1614601 1615681 1616658 1617828 1618921 1621068 1621081 1622161 1623241 1624308 1625295 1625388 1625401 1626373 1626482 1628628 1628644 1629725 1631770 1631886 1631888 1634028 1635128 1636088 1637268 1639430 1639448 1642563 1643750 1644848 1645800 1646990 1647008 1650230 1650248 1651197 1651311 1652276 1653488 1654551 1654568 1656728 1657671 1657808 1658871 1658888 1660910 1661048 1662126 1663193 1664148 1664286 1665353 1665366 1666305 1668593 1668606 1669673 1670766 1671701 1671833 1675085 1676018 1676165 1677233 1678175 1678325 1679393 1679405 1680485 1681411 1682633 1684803 1685873 1688043 1688965 1691122 1691283 1692353 1692363 1694523 1695438 1695593 1696515 1697753 1698673 1698843 1701003 1701910 1702073 1702083 1703163 1705313 1705322 1707473 1708382 1708562 1710540 1710713 1712697 1712882 1715042 1716122 1717014 1717193 1717202 1718273 1719362 1720248 1720442 1722602 1723486 1723611 1723613 1725642 1725768 1725774 1725842 1726847 1727993 1728878 1729003 1729016 1729082 1730162 1731035 1732113 1732242 1732322 1733393 1733402 1734270 1735349 1735561 1735562 1736554 1736633 1736642 1737505 1737721 1737722 1737723 1738582 1738741 1738803 1738804 1739792 1739823 1739882 1740962 1740965 1741949 1741983 1742047 1742897 1743122 1743128 1745054 1745187 1745290 1746306 1746363 1746373 1749369 1749549 1749593 1749614 1750448 1750629 1750683 1750696 1751662 1752605 1752789 1752833 1752860 1753869 1753923 1753941 1755843 1755977 1756083 1757107 1757184 1758233 1759213 1759323 1759346 1760156 1760405 1761233 1761370 1762587 1763392 1763586 1763645 1763669 1765805 1766628 1766824 1766885 1767844 1767965 1767990 1769072 1769861 1771205 1772019 1772158 1772273 1772313 1773096 1774475 1775253 1775462 1775525 1775556 1776328 1776331 1776605 1776636 1776638 1778485 1778630 1778753 1778799 1779845 1779881 1780643 1780925 1782042 1785180 1785245 1786185 1786364 1787115 1787445 1788473 1789269 1789270 1789565 1791579 1791658 1791767 1792502 1792805 1794953 1795010 1795734 1795737 1798253 1799131 1799218 1799285 1800046 1800414 1801125 1801288 1802513 1803277 1803657 1804739 1805698 1805765 1806509 1807760 1808668 1809062 1809743 1810142 1810143 1811225 1811899 1812177 1812233 1812245 1813150 1813386 1814056 1814228 1814468 1815471 1815485 1816495 1816629 1816631 1817289 1817712 1818365 1818540 1818542 1819791 1820955 1821599 1822037 1822853 1823755 1823756 1823931 1824200 1825191 1825205 1826083 1826213 1826987 1826988 1827444 1828239 1828373 1828431 1828529 1829143 1829511 1830605 1830690 1831470 1831775 1832375 1832692 1833625 1834529 1834530 1834704 1835932 1836005 1836102 1836855 1837183 1837761 1837764 1838151 1839009 1840309 1840993 1840994 1841329 1841510 1842591 1843150 1843311 1843316 1843567 1844569 1844647 1844754 1845305 1845460 1845647 1845836 1846382 1846789 1847614 1848536 1848689 1848887 1848949 1850165 1850690 1850693 1850841 1851249 1851764 1852990 1853287 1853412 1853921 1854349 1854493 1855144 1856076 1856447 1857737 1858233 1858368 1858374 1858687 1859748 1859903 1860380 1860387 1860517 1860767 1862067 1862534 1862664 1863150 1864231 1864683 1864686 1864819 1865148 1866248 1866395 1866397 1866837 1866962 1866965 1867901 1867907 1868408 1868560 1869116 1869407 1869644 1869647 1870047 1870050 1870184 1870187 1871808 1871811 1872191 1872201 1872339 1872893 1872896 1872897 1873257 1873259 1873266 1873786 1874484 1874489 1874493 1875058 1875408 1875409 1875413 1876635 1876640 1877048 1877228 1877231 1877551 1877556 1877559 1877560 1878047 1878316 1878319 1878540 1878543 1878546 1878551 1878554 1878563 1878566 1878568 1878574 1878580 1878584 1878592 1878596 1878599 1878603 1878607 1878611 1878614 1878619 1878622 1878623 1878626 1878628 1879186 1880482 1880678 1880681 1880686 1880689 1880692 1880697 1880942 1881567 1881570 1881571 1882007 1882016 1882801 1882808 1882811 1882816 1882825 1882828 1882832 1883446 1883528 1883737 1883744 1883746 1883752 1883757 1883761 1883767 1883770 1883780 1883788 1883793 1883802 1883808 1883817 1883820 1883825 1883829 1883834 1883837 1883840 1883844 1883848 1883852 1883855 1883857 1883860 1883863 1883864 1883867 1883872 1883878 1884151 1884161 1884583 1886297 1886306 1888447 1890008 1890585 1890586 1890589 1890596 1890599 1891616 1891629 1891630 1891638 1891643 1891651 1891657 1893166 1894328 1894826 1894828 1894834 1894839 1894845 1894848 1895381 1896975 1896981 1899125 1899128 1899646 1899728 1900200 1900202 1902357 1904514 1905589 1906181 1907744 1908286 1908368 1909421 1909898 1912052 1913686 1913768 1914206 1914766 1914848 1915283 1915899 1920365 1920374 1920382 1920384 1920387 1920391 1920397 1920398 1922326 1922379 1922522 1922564 1922826 1923406 1923488 1923591 1923594 1923648 1923897 1923898 1925648 1925746 1925816 1925819 1926646 1926728 1927778 1927898 1927989 1927990 1928206 1928209 1931048 1931132 1931237 1932098 1932325 1932330 1933126 1933208 1933289 1933594 1933597 1934206 1934288 1934496 1935745 1935748 1936519 1936670 1936673 1936816 1936821 1936824 1937446 1937528 1938578 1938838 1938844 1938849 1938856 1938858 1938872 1938877 1938888 1938901 1938906 1938914 1938925 1938930 1938938 1938945 1938957 1938965 1938971 1939109 1939120 1939131 1939141 1939154 1939165 1939170 1939175 1939183 1939188 1939195 1939200 1939210 1939606 1940686 1940836 1941375 1941766 1942896 1942928 1942992 1943418 1943926 1945718 1946229 1947166 1947248 1947889 1948246 1949376 1950061 1950488 1950545 1950949 1952236 1953646 1953728 1953782 1954407 1954415 1954863 1955856 1956337 1956585 1956591 1956887 1957025 1958485 1959047 1960268 1960930 1960935 1961255 1961288 1961349 1963101 1963367 1963512 1963514 1964528 1965268 1965277 1965575 1966688 1967687 1967838 1968529 1969615 1969630 1969847 1969928 1970003 1970333 1970926 1970930 1970932 1971089 1971091 1971195 1971216 1971233 1971255 1971265 1971275 1971281 1971289 1972055 1972386 1972388 1972400 1972405 1972471 1972482 1972896 1972917 1972937 1972985 1973026 1973044 1973061 1973261 1973264 1973266 1973269 1973270 1973273 1973285 1973295 1973300 1973308 1973319 1973340 1974570 1974575 1974583 1974588 1974604 1974615 1974623 1975249 1975295 1975328 1978489 1978568 1979569 1979648 1982888 1984969 1985015 1986049 1988209 1988255 1988288 1989289 1989335 1989367 1992607 1995815 1996927 1999055 2000133 2000167 2004453 2004485 2005533 2006645 2009885 2010932 2014172 2014205 2017412 2017445 2019570 2019605 2021222 2021225 2021233 2021238 2021243 2021250 2021259 2021269 2021281 2021286 2021289 2021302 2022810 2022845 2023347 2023350 2023354 2023360 2023363 2023373 2023376 2023468 2023470 2023476 2023478 2024371 2024385 2024390 2024395 2024403 2024404 2024409 2024419 2024566 2026049 2026509 2026512 2026520 2026525 2026729 2027165 2028646 2028649 2028654 2028663 2028664 2029978 2029982 2030777 2030786 2030792 2030797 2031484 2032146 2032149 2032919 2032926 2032932 2033230 2033605 2033971 2033976 2033979 2033992 2035395 2035399 2035400 2036102 2036105 2036107 2036115 2036118 2036122 2036844 2037566 2037964 2038253 2038257 2038648 2038651 2038653 2038656 2039002 2041478 2041482 2041900 2041903 2042986 2043320 2043364 2044397 2044699 2044701 2044707 2045159 2045524 2046553 2046556 2046604 2046823 2046828 2046836 2046839 2046842 2046845 2047321 2047325 2047331 2047631 2048764 2048970 2048975 2049498 2049502 2049504 2049787 2049845 2050592 2050593 2051121 2052761 2053021 2053023 2053088 2053273 2053278 2054344 2054345 2054924 2055173 2055178 2055252 2056012 2056013 2056016 2056249 2056333 2056336 2056497 2056500 2058181 2058186 2058187 2058192 2058196 2058409 2058499 2058502 2058505 2058508 2058510 2058516 2058517 2058520 2058522 2058648 2058651 2058652 2059277 2059282 2059440 2059441 2059445 2059448 2059452 2059455 2059458 2059459 2059461 2059464 2059465 2059466 2059470 2059491 2059492 2059605 2059606 2059609 2059613 2059617 2059618 2059621 2059624 2059626 2059629 2059715 2059721 2059724 2060574 2061453 2061594 2061595 2061596 2061635 2061637 2061646 2061649 2061663 2061667 2061795 2061858 2061860 2061863 2061866 2062668 2062669 2062670 2062737 2063617 2063625 2063631 2063635 2063957 2063960 2064012 2064015 2064825 2064826 2064828 2064902 2065052 2065057 2065061 2065065 2065070 2065075 2065078 2065080 2065084 2065087 2065091 2065813 2065905 2065985 2067978 2067980 2067993 2068001 2068002 2068007 2068010 2068016 2068063 2068146 2068149 2069140 2069142 2070183 2071265 2071269 2071272 2071276 2071279 2071284 2071287 2071290 2071296 2071298 2071299 2071300 2071397 2073465 2073470 2073475 2073560 ];
L = superpixels(X,10468);
BW = lazysnapping(X,L,foregroundInd,backgroundInd);

% Create masked image.
maskedImage = X;
maskedImage(~BW) = 0;
end

