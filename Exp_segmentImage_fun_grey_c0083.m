function [BW,maskedImage] = Exp_segmentImage_fun_grey_c0083(X)
%segmentImage Segment image using auto-generated code from imageSegmenter app
%  [BW,MASKEDIMAGE] = segmentImage(X) segments image X using auto-generated
%  code from the imageSegmenter app. The final segmentation is returned in
%  BW, and a masked image is returned in MASKEDIMAGE.

% Auto-generated by imageSegmenter app on 02-May-2022
%----------------------------------------------------


% Graph cut
foregroundInd = [574288 574291 574296 574299 575365 576464 581840 581866 586160 587266 588320 588346 588350 593720 594800 594837 595917 600200 601317 601322 607760 609968 609969 612080 614294 617480 619640 619699 622946 625045 628349 630445 631525 632669 632672 635919 635920 636928 640173 640240 642405 644493 648891 648892 649893 653139 653212 655377 655379 660779 663939 665019 666179 666184 667264 670419 672579 672669 675899 675902 675909 677979 678062 682299 683471 684459 688871 688873 689859 691038 695259 702819 704000 708219 710373 710486 712646 714693 717933 718050 722253 722370 725613 725617 727653 729813 731017 736422 741693 743848 745064 754648 756808 760189 761272 765448 769768 773157 774237 779488 783808 792445 793685 796925 800005 808810 810805 811885 815295 822855 823760 827000 828255 831495 835640 835816 840142 844280 847517 850947 851837 856347 856350 859397 862830 865873 868234 871273 877753 877957 879907 883357 886387 887682 891786 891787 892002 895242 898266 903889 906906 910370 918780 919010 927420 927650 930655 934136 934975 937376 939295 942532 946020 953580 954412 957647 958980 966540 969527 970860 977083 978160 982480 989224 993280 998675 998944 1001915 1012710 1018384 1021624 1025668 1033222 1041064 1041862 1045384 1053739 1054024 1064539 1065615 1065904 1072095 1077490 1081024 1083967 1085344 1087504 1089362 1091522 1096917 1101235 1104470 1108785 1111260 1113103 1113105 1116338 1122818 1123140 1124975 1125294 1131451 1135014 1136846 1137923 1143323 1145478 1149054 1150874 1154449 1160591 1163089 1164905 1168487 1170305 1173544 1176047 1185418 1185762 1187573 1191162 1197290 1198716 1202684 1208084 1210239 1210596 1215637 1221032 1223192 1223553 1227507 1229666 1233268 1235060 1240828 1241537 1244773 1250168 1254485 1254866 1256640 1256645 1263117 1263506 1266357 1269593 1274989 1276069 1277541 1282546 1286181 1287942 1290096 1295496 1298061 1300216 1300895 1304129 1312096 1314256 1317081 1319654 1320729 1320734 1322476 1328284 1330033 1331523 1331524 1334757 1335429 1342985 1347708 1348382 1349864 1351616 1355261 1355931 1355932 1359171 1364565 1364975 1369960 1370370 1371037 1376437 1379001 1381832 1385069 1388716 1389389 1394114 1394785 1396938 1396941 1398428 1402338 1405568 1405574 1408148 1409227 1409888 1415283 1415287 1418517 1418523 1423910 1423916 1423917 1424336 1427150 1429733 1430389 1434048 1435782 1437938 1443338 1443762 1448735 1449160 1449815 1451310 1451315 1456290 1456704 1456709 1456710 1459930 1459933 1459938 1459941 1461690 1464925 1465321 1469626 1469633 1472862 1474645 1475723 1478248 1478259 1485794 1489763 1491191 1493346 1498398 1499823 1503798 1505215 1505218 1508452 1509193 1509198 1511353 1511687 1515673 1517083 1517087 1520320 1525391 1525692 1525697 1525702 1525708 1525710 1530786 1531066 1531073 1531077 1531079 1531090 1533213 1533218 1539422 1539688 1542920 1544822 1545902 1545905 1545909 1545911 1545916 1545923 1545927 1545931 1545934 1545939 1545946 1545952 1545963 1547235 1551363 1553706 1556770 1556774 1558934 1563418 1564338 1568810 1569738 1569741 1569745 1569752 1569754 1569763 1569765 1569771 1569774 1569776 1569785 1569789 1569793 1569800 1569805 1569812 1569818 1569826 1570909 1570913 1570961 1570965 1576358 1577400 1577409 1579571 1579575 1579580 1579581 1579586 1579588 1579593 ];
backgroundInd = [107552 107559 107581 107590 107602 107607 107619 107624 107631 111827 111848 111851 111859 111960 111965 116120 116127 116140 119533 119552 120412 120419 120422 120433 120640 120641 120655 126865 126872 126880 127143 127153 128981 128992 128999 129003 129012 129019 129324 133797 133799 133806 133810 133819 133823 133828 134292 134314 134321 134335 134340 134349 134359 134365 134370 134378 134725 134736 134746 139228 139680 139682 140150 140166 140173 140184 140186 140194 140212 140217 140224 140229 140246 140250 140257 140266 141390 145710 145714 150441 150449 152197 152202 153282 153289 155830 159061 159774 164094 164100 164456 166260 169854 170934 171661 174901 177409 178141 180608 186049 191447 192527 194338 206567 213757 222767 235727 242897 243287 269201 269883 278519 278921 282161 295784 300521 307650 311321 317364 324903 324908 326441 332447 335081 338922 340481 343721 346961 348623 359921 360488 362646 366958 373961 376667 379361 391241 391773 398801 399324 406882 415509 423067 425801 433354 438754 454954 467912 481349 481952 484112 492752 498152 502472 504021 509416 511112 532707 534246 541347 542886 557995 573747 577433 580668 589308 601185 611541 615861 626013 641131 654737 661645 672438 679996 689290 695111 696191 701586 706984 715619 732890 745844 748684 749083 753403 762037 771752 777842 786872 795115 798748 811315 823582 846953 847334 863148 875404 879342 893758 906713 918220 922911 941260 946658 948460 960689 970404 984429 998462 1005693 1008176 1021123 1031920 1044865 1044869 1056741 1067534 1084807 1089127 1092086 1107479 1120432 1124482 1139865 1154979 1174410 1176565 1181716 1198163 1203558 1221909 1227305 1234861 1235713 1241336 1252987 1258607 1262709 1268323 1280197 1282354 1307179 1328775 1334172 1341726 1347126 1394628 1400023 1420537 1432417 1448613 1465888 1471288 1491802 1497202 1520962 1530682 1555518 1567398 1585204 1585206 1585210 1585215 1585223 1585230 1585239 1585243 1587346 1587353 1587411 1587422 1591158 1592824 1593820 1599213 1599313 1600291 1600393 1605680 1605685 1605793 1606278 1608913 1608918 1609038 1613229 1620786 1622004 1628342 1630653 1631114 1634812 1636964 1639297 1642354 1645592 1652061 1652714 1657458 1659610 1660906 1663068 1666077 1668914 1671475 1672788 1674704 1679015 1679017 1687919 1688723 1695201 1698436 1701314 1702750 1705203 1709523 1711034 1711387 1725423 1730474 1740194 1744518 1752408 1753161 1760726 1766126 1768933 1775854 1778653 1784798 1788820 1791613 1797466 1810434 1815373 1819077 1819336 1821848 1826648 1826894 1834211 1834803 1841783 1842008 1843443 1850426 1850644 1850999 1853670 1858194 1866640 1866826 1867195 1868805 1874212 1874382 1878692 1880701 1881226 1881786 1881927 1882301 1888272 1888393 1888398 1888778 1890439 1890442 1890542 1890545 1894166 1896931 1896956 1897406 1902357 1902370 1902379 1902392 1902400 1902403 1902726 1902728 1902746 1902757 1902760 1902771 1902778 1902783 1909180 1909190 1910244 1910248 1910255 1917787 1917794 1922102 1932885 1937654 1937665 1941520 1944139 1945832 1949545 1949546 1949550 1956619 1957119 1957127 1961452 1962017 1965248 1969550 1971191 1976594 1977104 1979839 1979847 1985251 1985253 1985263 1985266 1985271 1985278 1985283 1985290 1985299 1985302 1985316 1985323 1985343 1985358 1985366 1985375 1985389 1985428 1985435 1985443 1985456 1985463 1985476 1985485 1985734 1987655 1987665 1987669 1987687 1987698 1987699 1987717 1987731 1987875 1987883 1993143 1993153 1993164 1993171 1993186 1993193 1993204 1993211 1993216 1993223 1993236 1993254 1993263 1993270 ];
L = superpixels(X,10468);
BW = lazysnapping(X,L,foregroundInd,backgroundInd);

% Create masked image.
maskedImage = X;
maskedImage(~BW) = 0;
end

