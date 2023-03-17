function [BW,maskedImage] = Exp_segmentImage_fun_RGB(RGB)
%segmentImage Segment image using auto-generated code from imageSegmenter app
%  [BW,MASKEDIMAGE] = segmentImage(RGB) segments image RGB using
%  auto-generated code from the imageSegmenter app. The final segmentation
%  is returned in BW, and a masked image is returned in MASKEDIMAGE.

% Auto-generated by imageSegmenter app on 25-Feb-2022
%----------------------------------------------------


% Convert RGB image into L*a*b* color space.
X = rgb2lab(RGB);

% Graph cut
foregroundInd = [570983 576387 578547 578548 578553 578555 578559 583959 583964 583966 589371 592616 592626 595840 598029 598031 600200 606691 612093 617443 620746 626151 629397 634804 640208 643375 643458 648861 648865 652105 655350 655352 660757 663999 664003 668328 671456 675897 675902 676983 683469 688874 688875 691044 696450 698613 704017 704021 706022 709421 712502 712666 714662 718066 718986 721308 724386 724548 726546 729953 737346 738595 744906 748323 753726 754626 761286 766689 772089 772094 774259 779661 781826 783786 787026 788311 792631 792633 800193 804518 807759 808626 812079 819645 820730 822666 823970 823973 834773 839940 840178 847744 852064 853980 862615 871511 875575 875831 879895 892031 898252 903652 903911 911471 915791 919031 923086 927401 927671 934961 941441 941711 946839 948994 951431 954394 966266 970871 974897 978431 982452 984911 993244 993551 998639 998951 1006190 1009751 1012662 1021297 1021631 1025608 1036402 1037477 1037831 1043955 1045386 1049355 1052590 1056905 1065543 1068066 1073094 1076706 1080649 1085346 1089283 1097226 1101161 1104396 1104786 1111266 1116270 1120986 1124907 1128142 1132866 1140022 1143260 1143664 1154055 1157291 1158778 1164848 1167418 1170242 1179953 1181453 1184693 1189013 1192907 1193333 1196573 1200890 1202622 1208020 1210605 1215999 1216651 1226365 1233277 1238672 1246870 1247308 1251183 1260256 1263058 1265212 1269974 1272129 1286160 1286803 1288961 1293715 1302349 1307741 1309469 1317458 1322847 1324574 1325006 1327810 1330406 1335800 1337955 1340767 1342916 1344427 1353064 1354792 1359533 1366658 1368170 1373127 1376805 1378524 1383276 1388236 1391911 1394703 1397303 1402261 1407018 1408729 1411966 1412409 1417363 1419517 1421041 1424910 1424914 1430308 1430310 1431830 1432458 1437220 1437856 1442619 1443246 1443251 1444773 1445404 1445406 1450800 1451248 1452957 1457273 1459873 1465266 1466986 1468491 1473879 1478864 1481019 1483586 1483597 1490057 1492899 1493292 1496139 1496522 1496530 1499379 1502624 1504070 1510538 1515929 1518839 1519157 1519162 1525326 1525624 1525629 1530743 1531003 1531008 1531016 1534003 1534007 1534211 1534215 1534223 1534227 1534236 1539421 1539428 1539439 1539447 1539459 1539466 1539479 1539487 1539492 1539499 1539505 1539512 1539520 1539538 1539542 1539550 1539553 1539564 1539571 1539582 1539589 1539604 ];
backgroundInd = [451023 451027 451032 451034 451038 451043 451051 451052 451062 451064 451069 457495 457500 457556 457562 457567 457573 462890 462979 462983 462986 462990 462994 462997 463001 463003 465044 465048 465163 465167 465172 465175 465178 470583 470591 470593 470597 473681 473842 473851 479262 484667 484670 486636 486835 488791 498730 500667 500895 500904 506304 509549 509554 510382 512537 514956 517937 520361 520362 522522 526574 527928 533336 536294 536588 539534 541990 542774 547398 548174 550642 551414 556814 557131 559302 564374 564704 569774 571197 577337 580919 585982 586324 590302 595702 598942 600373 604347 605427 607937 616230 618390 623065 629194 631354 631711 636759 637113 639999 643241 646845 652961 664841 665213 670616 673487 677807 678181 684666 687529 689689 692226 694009 698329 701574 701949 706974 712374 714539 716000 719939 720320 727499 730043 732899 736139 744088 745861 748021 757741 766771 767461 768541 776495 779341 779735 787981 792700 793381 799861 807421 812140 814981 816061 820784 831181 835909 843061 844141 844549 851701 852109 858179 865739 866149 871139 875874 879779 880194 883434 887754 891654 892074 899634 901374 903528 906114 906768 911514 917568 922314 922965 927281 927714 930521 934833 935274 941313 946074 947154 951020 965055 966128 966594 969834 978474 981240 982313 982794 987114 994674 998498 998994 1001154 1006056 1006554 1009794 1012531 1013034 1017924 1021158 1021674 1030314 1033034 1033554 1036269 1037874 1041665 1043822 1044354 1049218 1053531 1056763 1065954 1068640 1072950 1073509 1077269 1081069 1085389 1089135 1089709 1092949 1097269 1103749 1110717 1111309 1112872 1116709 1119944 1120425 1125344 1127978 1127979 1131824 1132288 1132293 1134445 1137224 1139384 1139837 1139840 1148468 1149104 1151264 1153866 1156664 1160335 1162064 1162492 1163144 1167887 1173279 1175434 1182584 1184744 1185147 1189064 1194462 1197702 1198097 1200252 1205262 1207806 1215359 1216062 1218222 1220751 1220754 1222542 1227942 1228306 1230102 1234779 1238013 1243062 1248457 1249885 1252040 1258177 1258517 1266071 1273626 1274706 1276537 1278697 1283017 1286257 1286579 1290577 1295215 1295971 1300610 1309250 1311088 1317885 1318648 1328363 1331917 1332683 1340552 1341321 1349955 1353508 1355350 1363223 1367543 1377259 1378021 1382656 1385576 1387736 1390216 1394214 1397449 1400689 1403170 1411810 1413642 1413968 1421197 1422274 1424768 1426928 1440963 1444940 1450334 1460049 1460398 1465446 1467955 1483071 1483793 1492428 1499271 1504298 1511144 1514013 1520484 1525184 1526264 1528039 1531664 1532350 1537064 1537748 1543143 1544624 1545297 1550691 1553269 1557159 1558236 1558669 1565780 1568395 1569009 1573798 1574406 1574878 1579793 1580278 1582438 1582442 1587842 1593248 1593249 1593254 1593256 1593260 1593811 1596505 1596507 1596512 1596517 1596518 1596523 1597040 1608907 1611058 1617531 1620758 1631544 1634776 1640849 1641237 1641248 1646624 1648430 1648762 1648769 1648774 1648777 1654144 1654149 1657368 1658171 1658173 1663592 1663611 1663632 1663815 1663823 1663828 1663837 1669045 1669058 1669059 1669066 1669070 1669079 1669083 1669087 1669098 1669109 1669125 1669153 1669171 1669175 1669191 1669193 ];
L = superpixels(X,10468,'IsInputLab',true);

% Convert L*a*b* range to [0 1]
scaledX = prepLab(X);
BW = lazysnapping(scaledX,L,foregroundInd,backgroundInd);

% Create masked image.
maskedImage = RGB;
maskedImage(repmat(~BW,[1 1 3])) = 0;
end

function out = prepLab(in)

% Convert L*a*b* image to range [0,1]
out = in;
out(:,:,1)   = in(:,:,1) / 100;  % L range is [0 100].
out(:,:,2:3) = (in(:,:,2:3) + 100) / 200;  % a* and b* range is [-100,100].

end
