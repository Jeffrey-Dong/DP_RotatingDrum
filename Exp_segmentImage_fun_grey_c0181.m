function [BW,maskedImage] = segmentImage(X)
%segmentImage Segment image using auto-generated code from imageSegmenter app
%  [BW,MASKEDIMAGE] = segmentImage(X) segments image X using auto-generated
%  code from the imageSegmenter app. The final segmentation is returned in
%  BW, and a masked image is returned in MASKEDIMAGE.

% Auto-generated by imageSegmenter app on 27-Jul-2022
%----------------------------------------------------


% Graph cut
foregroundInd = [553870 553871 553874 553876 553879 554946 554948 557104 558184 558202 559287 560343 563581 563610 564694 566820 569022 571138 571184 573297 573345 575508 576535 576590 578753 579775 579834 581933 583076 585172 585237 585239 586252 586319 588483 590570 591650 591725 595970 597049 597128 600289 600370 602446 604692 605686 607935 608926 611175 616485 616577 618741 619723 620803 621983 625121 625225 627386 628361 630626 632680 633760 634948 637109 638080 643478 643591 645752 646716 648876 648994 650075 653195 654395 656435 658717 661833 661960 663992 664121 667232 667363 670603 672764 673709 676004 678027 680326 681267 681407 682346 684647 686809 689904 690050 692063 694370 695452 696381 698692 699619 700853 702859 704095 705018 707337 709499 713820 714734 717062 717971 719051 719224 722465 723370 725707 726608 727687 730029 732005 732189 734164 735430 737402 740640 742799 744074 746037 746234 747116 749477 751436 755751 755959 758990 760068 760280 761360 764387 764602 766763 767625 768704 771083 774102 774325 776485 777340 779726 780579 780808 782737 782968 784048 785976 787289 791373 792691 794611 795690 795934 798928 799176 801087 803496 804325 804576 807564 807817 810804 811057 814042 814297 815121 817539 818619 819439 821859 822676 824019 825915 827259 828073 828340 832660 835900 836710 837788 840220 841027 841302 844542 845345 846424 846703 850742 851023 854263 855060 855345 856139 858585 859376 860745 863695 863986 864773 867228 868011 871548 875568 875870 878806 879110 880192 880965 883432 884205 885283 885592 889914 890994 893919 894234 894997 896394 899635 900715 901474 902875 903955 904712 906117 910110 910437 913677 915508 916917 919827 920157 923397 925223 927717 929877 932778 933117 936017 936357 938518 941758 942838 943573 944998 945733 949318 950051 952558 953638 956528 959038 959766 960118 963006 964438 968758 969838 970561 973078 973799 975238 976318 977398 979558 980278 980638 982798 983516 983878 984595 986038 987118 991438 993231 994678 995390 997918 998628 1002238 1002943 1003318 1004399 1006182 1008340 1008719 1010879 1011959 1012659 1015199 1015895 1019134 1019519 1020599 1021292 1023839 1024531 1025609 1028159 1031008 1031399 1032085 1035324 1035719 1036799 1037482 1040039 1040721 1043279 1043959 1045439 1046118 1049759 1050434 1050839 1051919 1053672 1054079 1055159 1055831 1058399 1059068 1060559 1061639 1062307 1064467 1068119 1068785 1069864 1071359 1074599 1076340 1078919 1079578 1080657 1082159 1084977 1088639 1090799 1091453 1094039 1095119 1095772 1098358 1099009 1100518 1101598 1104838 1105485 1109158 1109802 1113040 1113476 1114556 1115199 1117796 1118437 1118876 1121675 1123196 1123834 1124276 1125356 1128149 1128596 1129676 1132468 1132916 1135706 1136156 1136785 1140476 1141103 1141556 1144341 1144796 1146500 1146956 1149738 1151276 1151897 1154514 1155132 1157754 1159451 1159914 1160994 1162074 1162687 1163766 1166394 1169164 1172401 1172874 1173480 1176114 1176718 1178877 1179353 1182112 1184753 1187510 1189073 1190153 1190747 1191826 1194470 1196142 1198790 1199381 1201536 1202030 1204190 1206934 1208510 1209590 1210171 1210670 1212830 1213408 1213910 1215566 1216070 1218803 1219310 1222550 1223120 1225277 1226868 1228514 1229028 1230672 1232268 1233909 1235508 1238227 1238748 1239304 1243067 1243622 1244147 1247941 1248465 1249018 1251705 1252255 1254413 1256025 1257104 1257652 1258730 1260344 1263584 1266824 1267367 1267902 1269526 1271142 1272764 1273302 1273843 1276541 1279779 1282478 1285714 1287337 1288417 1291656 1292190 1294896 1295429 1298667 1299214 1305691 1306220 1307851 1310539 1311090 1312170 1315410 1317568 1319175 1320251 1320807 1325125 1325649 1328888 1329443 1329965 1330522 1334284 1335918 1337522 1338601 1339157 1343999 1344554 1345076 1345633 1348871 1353189 1353712 1355868 1356428 1358586 1359106 1361825 1362903 1364504 1366660 1367222 1370462 1372620 1374217 1375296 1375859 1376937 1378534 1381255 1385573 1386090 1389327 1389889 1390968 1393644 1394206 1394723 1396366 1397961 1399605 1400683 1403358 1404436 1405002 1408753 1409317 1410396 1411988 1413633 1415227 1417385 1417949 1420622 1421187 1426020 1426586 1427097 1428742 1431415 1431978 1433570 1435216 1436809 1437373 1441691 1443847 1444363 1447599 1448166 1449244 1451918 1452483 1452993 1454640 1457310 1457877 1458389 1460034 1463785 1464862 1468667 1469180 1471905 1472415 1474572 1476222 1477301 1477810 1480539 1483778 1485363 1486442 1487016 1489679 1490254 1491836 1492411 1495649 1497229 1499388 1499968 1503704 1504286 1505861 1508016 1508601 1509680 1511253 1512332 1512918 1515077 1516648 1518315 1519394 1519885 1522043 1523712 1525869 1526359 1528515 1530184 1531754 1533421 1535577 1538228 1538814 1541464 1542052 1544700 1545289 1545779 1547445 1550095 1550683 1551172 1552840 1553331 1556078 1556567 1558235 1558723 1561960 1562553 1563038 1564709 1565787 1567355 1568431 1570104 1570588 1571180 1572744 1573819 1573822 1574898 1575499 1577055 1577655 1578133 1581370 1581971 1583049 1584606 1586764 1587366 1587841 1588443 1589997 1592761 1594314 1594918 1595390 1598156 1600788 1601391 1601864 1604628 1606183 1608339 1608944 1610023 1611576 1613260 1615416 1615894 1616970 1620814 1621890 1622368 1625606 1626207 1626209 1627763 1629446 1630522 1631001 1632679 1636398 1636994 1638073 1638556 1641310 1642875 1643955 1644546 1646112 1646703 1647192 1649941 1650431 1652096 1652591 1654253 1655829 1657491 1660149 1660728 1661228 1662885 1663961 1665548 1666119 1667196 1668786 1669352 1670431 1672024 1672588 1674184 1674744 1675820 1677423 1677976 1681211 1681741 1682289 1683901 1684446 1685525 1687139 1687681 1688759 1690379 1690915 1690916 1691458 1693072 1694144 1694148 1695778 1696301 1697376 1697379 1700098 1700611 1700613 1701178 1702768 1703843 1703846 1704418 1704921 1706578 1707078 1709232 1709235 1709818 1710307 1712464 1713058 1714615 1714618 1715694 1716298 1720618 1721082 1721086 1724938 1725395 1725397 1726017 1726019 1726471 1727099 1727100 1728178 1728626 1729258 1729260 1729261 1730339 1730341 1730782 1732499 1732503 1732504 1732933 1732937 1733579 1733584 1734005 1734007 1734008 1734011 1734666 1736156 1736159 1736161 1736827 1736828 1737231 1737234 1738988 1738990 1740464 1740467 1741150 1741151 1741535 1741539 1741541 1742234 1743690 1743692 1744395 1744397 1744398 1744399 1744401 1745484 1745485 1745489 1745493 1745844 1745847 1746917 1746920 1746923 1747654 1747657 1747658 1747661 1747663 1748746 1748751 1748754 1748755 1749073 1749074 1750149 1750918 1750920 1750923 1750926 1752007 1752011 1752015 1752305 1752306 1754178 1754179 1754182 1754184 1754188 1754189 1754192 1754457 1754461 1755273 1755280 1755281 1755285 1755287 1755531 1755534 1757455 1757458 1757459 1757462 1757689 1759843 1760710 1760711 1761997 1762000 1762875 1762878 1763959 1763963 1763965 1763969 1763971 1763974 1764154 1765055 1765060 1765063 1765064 1765067 1765072 1765076 1765078 1765081 1765225 1765226 1767245 1767250 1767253 1767257 1767374 1767375 1767380 1769419 1769422 1769425 1769430 1769433 1769436 1769514 1769517 1769519 1769522 1769525 1769527 1769531 1770519 1770522 1770528 1770531 1770535 1770536 1770539 1770541 1770544 1770547 1770548 1770553 1770556 1770559 1770560 1770563 1770566 1770568 1770572 1770573 1770579 1770582 1770589 1770591 ];
backgroundInd = [2276 2280 2281 2284 2286 2289 2290 2293 2296 2299 2302 2303 2306 2308 2311 2312 2315 2319 2321 2324 2327 2329 2332 2333 2341 2344 2346 2349 2352 2355 2359 2361 2364 2367 2368 2372 2374 2377 2382 4349 4351 4352 4354 4355 4357 4360 4361 4364 4368 4370 4373 4374 4377 4378 4382 4386 4387 4390 4394 4397 4400 4401 4405 4408 4410 4413 4416 4419 4423 4425 4428 4429 4432 4545 4547 4551 4554 4556 4559 4562 4563 4566 4569 5427 5652 5657 5661 5667 5670 5673 5675 7584 7838 7841 7842 7850 7851 7855 8496 8498 8501 8505 8506 8509 8512 8514 8517 8520 8524 8526 8529 8532 8533 8537 8541 8542 8545 8546 8550 8553 8554 8557 8561 8564 8566 8569 8570 8662 9535 9538 9539 9542 9545 9548 9551 9552 9555 9561 9564 9566 9569 9653 9655 9659 10020 10021 10025 10821 11112 11113 11116 11120 11126 11129 11654 11657 11660 11662 11665 11673 11674 11677 11678 11682 11685 11686 11825 11828 11900 12978 13292 13294 13300 13307 13308 13783 13786 13792 13796 13801 13804 13808 13991 13994 15079 15137 15472 15475 15479 15480 15483 15489 15492 15493 15497 15501 15920 15923 15925 15928 15931 15932 15936 15938 16215 16584 16585 16588 16591 16593 16597 16598 16605 16606 16609 17240 17243 18066 18067 18070 18076 18324 19856 19859 19862 19867 19869 19875 19880 19886 19899 19905 19912 19914 19921 19927 19931 19934 19939 19942 19945 19952 19955 19956 20196 20206 20209 20210 20217 20219 20535 21568 22123 22126 22132 22135 22141 22142 22149 22155 22158 22162 22167 22171 22174 22178 22179 22183 22186 22187 22190 22197 22200 22201 22204 22209 22213 22218 22220 22225 22228 22235 22240 22243 22249 22257 22262 22266 22272 22279 22280 22288 22289 22296 22297 22304 22306 22309 22315 22318 22321 22327 22330 22332 22335 22338 22339 22347 22348 22355 22651 22652 24854 26977 28058 28092 29172 31301 33463 33489 34543 36704 36729 37785 41049 42107 44289 46429 49669 50749 50750 50769 52910 54008 56150 57230 58328 61550 61568 63713 66953 68048 69128 70193 72368 73434 74528 77754 78848 79928 83156 83168 84248 85316 86396 87488 88556 91796 91809 92876 95036 95049 98289 101514 101529 104754 104769 108009 108190 108208 108220 108234 109092 110412 111234 112332 112587 113049 113064 113074 113092 113100 113116 113122 113132 113143 113155 113171 113186 113569 115169 115187 115197 115358 115363 115374 115380 115838 116634 117732 118383 118394 118625 118631 118635 118640 118648 118658 118665 118670 118812 119096 119758 119766 120033 120525 121268 121932 121935 122052 123114 123443 123454 123750 124101 125292 125415 126263 126705 127346 127349 127452 128876 128888 129131 129511 130594 130597 130692 131060 131754 132956 133230 133242 133932 134517 134923 135012 135417 137087 137154 139250 139332 140974 141492 141584 141908 141919 143577 144091 144658 144732 146257 146268 146275 146283 146294 146299 146310 146315 146324 146331 146339 146352 146874 147972 148982 149132 150063 151212 152272 153372 155512 157692 157759 158706 160932 161946 162009 163072 164233 165249 166312 168550 169552 170588 171729 174952 174969 175026 179228 179272 181432 182590 183609 183666 183670 184672 185769 185830 189009 189070 190028 190072 190148 190150 192249 193329 193390 194390 195550 196508 197649 198710 199870 201950 205270 206289 206350 207350 207369 208388 208430 210609 211750 212750 213788 218169 219230 219312 221409 222470 225668 225710 226481 226484 226489 226497 226505 226515 226527 226532 226538 226544 226748 226809 226872 228630 228632 228636 228709 228718 228723 230049 230783 230888 230895 230902 231110 232935 232937 233064 233068 233268 233352 234153 234308 237398 237400 237588 238325 238481 239769 239832 240478 240643 240788 240849 243884 244966 244967 245148 246249 246314 247126 247270 248027 248350 248388 249489 250548 251590 252334 253787 253874 254678 254867 257049 259150 259187 259888 260289 260354 261369 262390 263507 264609 266836 267790 268708 268905 269985 270675 271030 271089 273225 274329 275182 278226 278649 278716 279705 280748 282736 284025 284116 284700 286146 286209 288372 289423 292663 294785 295932 295996 296765 298063 298092 299804 301265 303240 303423 303463 304572 305198 306703 307878 309511 309903 311052 311872 312103 313183 317065 317532 319621 319758 320503 320772 321823 322903 324613 326978 327252 328398 328922 330211 330421 331543 332652 333237 334523 336972 337758 337979 338023 338625 339103 339132 341360 342071 342939 344226 345302 345613 347247 347456 347743 347773 348779 349398 349610 350983 351013 351767 352160 354178 355870 357158 357463 358014 358019 359312 359653 360163 360703 361813 362179 362197 362204 362221 362239 362253 362263 362275 362285 362294 362305 364327 365780 367132 367538 367549 368263 368293 368362 369343 370453 371169 372583 372923 374773 376151 376560 376850 377983 380175 380459 380872 381223 382303 383415 385186 385486 385543 387703 387802 388007 389895 390159 390579 390943 392023 393135 394124 395295 396343 398535 398787 400602 401373 402451 402922 404983 405015 407079 407143 407980 407983 407988 407992 407996 408001 408002 408005 408156 408255 410133 410169 410315 411253 411495 412286 412814 413364 413414 413553 413623 415710 416483 416599 416662 416863 417747 417867 417975 418755 419287 419912 420016 420018 421989 422077 422164 422170 422295 422959 423343 423375 424146 424241 424246 424315 424318 425330 425336 425380 425386 425389 425503 426302 426686 426839 427501 427508 427513 427520 427527 427532 428775 430077 430905 431696 432015 434934 435225 435912 436335 437631 440327 440655 441705 441735 442620 442625 442629 442630 442633 442637 442638 442643 442644 442647 442650 442651 443027 444645 444774 444777 444814 444945 445847 445849 445851 445898 445901 447784 448004 448005 448065 449080 449081 449082 449295 449366 449505 450039 450227 450230 450345 451238 451240 451311 452195 452315 452317 453615 454554 455553 455554 455638 456631 457594 457707 457709 458672 458987 459015 459221 459866 459959 460946 461044 461817 461911 462227 463104 463207 463307 464183 464289 464415 464619 465149 465262 466547 467422 468500 468610 469694 470659 470867 471097 471737 471739 472048 474135 474773 474977 475095 475943 476056 477023 477136 477348 478339 478428 479294 479296 479535 480501 481668 482534 482975 483613 483743 484578 484824 484908 485771 486015 487068 488066 489011 489453 490089 490308 490962 491054 492388 494408 494655 495808 496566 496711 497868 498873 499805 500028 500883 501330 501755 502113 503194 503268 504006 504123 504348 505355 505455 506282 506730 507588 508597 509520 509677 510828 510855 511839 513721 514070 514918 516161 517242 518157 518415 518610 519237 519404 520197 521395 521630 522261 522475 522645 522710 523554 523727 525282 525294 525304 525307 525887 525950 526051 526968 527250 527431 527470 527874 528050 529190 529215 529585 529635 530211 530716 530993 531111 531293 531971 532191 532430 533959 534534 535615 535670 535890 536053 536511 537855 538670 538858 538910 540365 540708 541018 541019 541442 541910 542099 542150 542606 543230 544261 544530 545148 545342 545491 545755 545847 546422 546495 547079 547550 547912 548388 548584 549665 549930 550068 550423 550548 550790 551145 551628 551825 551870 553410 553986 553988 554384 555135 555947 556151 556190 557232 557617 558107 558314 558570 559430 559774 560139 560478 560535 560853 561108 561347 561559 561561 563722 563725 564585 564807 564808 564830 565050 565167 565665 566970 566971 566990 567015 567324 567452 568052 568055 568402 568905 570218 570221 570230 571065 571303 571306 571307 571308 571310 573095 573495 573690 573798 573932 575139 575385 576814 577035 578113 578623 580170 580272 582135 582943 584732 584970 585570 586183 586750 587263 587535 588343 588810 588909 593742 594015 594822 595290 595385 595647 596612 598060 599005 600219 600495 601770 601862 603459 604537 605250 607056 607776 608055 608600 609330 610880 613173 613738 614253 614535 617492 618205 618313 619050 621015 621810 622378 622889 626129 627207 627778 628029 628575 629848 632098 632475 632606 633684 634059 635480 636135 636328 636922 637742 640162 641535 642321 643980 645559 648208 648433 648798 649095 649384 650694 650830 650955 654908 655275 656655 658514 660086 660672 662055 663911 664303 664727 664989 666777 667455 668104 673046 673628 674706 675015 675408 677945 678756 680103 681060 682575 682956 683342 684921 684969 686985 687129 689055 689817 691215 691451 692787 694135 694817 697931 698453 698775 700032 700039 700198 701172 701572 702192 702254 703275 703278 703850 704175 705255 705435 705495 705734 707657 707660 708168 708676 708741 708822 710838 711988 711993 712000 712012 712026 712044 712486 712998 713565 713895 714975 716238 716239 716240 717766 717883 718308 718400 722535 723002 723278 723615 724881 727597 731361 731799 731915 732255 735152 737197 737310 737655 740003 740270 740549 744865 745311 745405 746297 748104 749073 749181 752967 753498 754302 755657 758177 758895 760529 763105 763211 764657 766450 767014 768606 772821 772925 773297 774002 774811 775656 776631 778320 779977 781059 781457 781557 783446 784097 784300 785875 786954 788622 790577 791271 792942 793429 795486 796668 797059 797479 797746 799426 801963 802064 802459 802551 802667 803750 805302 806381 808071 810435 811779 812393 814339 814913 815635 817175 817798 819739 820309 820414 822120 822572 823387 825361 825810 826443 826783 829459 829553 829684 829865 833366 834340 834861 836339 836604 837251 838332 839735 841341 842815 843080 844818 846319 848061 848369 849981 850369 851445 851607 852795 853313 854954 856461 856711 857002 857919 858192 860076 861861 863313 864667 865470 866442 866715 867261 867793 870769 873024 873304 874103 874383 876180 876426 876981 877338 878062 879233 879492 879495 881592 881937 882680 882899 883462 883804 884096 885058 885931 885940 885945 885953 885959 892102 893690 893812 894353 897049 898128 898582 900165 901822 903403 905682 907222 907841 910956 912159 913793 914782 916476 918022 919592 919714 920793 923422 924030 926066 926189 930507 930982 933619 933745 934222 935904 936473 939142 939622 942253 942862 943461 946699 948858 949809 950422 952673 953176 953662 956282 956413 957492 959062 961810 962757 965049 966127 968782 969862 970446 972111 972473 972606 973551 975844 977422 978948 979082 981241 982822 984479 988311 988663 988798 990821 993116 994194 994702 997432 998378 999457 1001750 1002262 1002829 1006069 1008091 1009822 1011467 1014231 1014566 1014706 1015222 1015784 1019023 1021181 1023200 1023862 1024419 1029815 1030342 1030428 1030754 1033053 1034132 1035742 1038450 1041687 1042222 1042766 1043701 1048164 1049868 1050862 1051255 1051400 1055716 1056262 1058954 1060033 1060968 1062191 1062742 1065430 1067442 1069306 1070828 1071904 1074622 1075143 1077302 1080391 1080539 1082182 1083776 1087014 1088743 1089023 1089740 1090253 1091331 1095140 1095497 1095649 1098888 1099966 1104129 1104285 1107020 1108603 1111684 1112420 1112500 1115076 1118158 1118313 1120472 1121060 1123710 1124628 1127539 1128027 1130185 1130860 1132182 1133423 1134019 1136661 1138820 1142056 1145899 1146209 1146375 1147058 1148533 1151299 1151771 1153759 1155010 1157168 1158859 1160407 1163643 1166419 1166498 1166881 1167787 1170119 1171819 1174438 1175514 1177219 1178574 1178753 1180911 1183699 1184148 1188465 1190255 1190624 1192339 1192601 1193861 1196019 1198817 1199075 1199257 1203574 1205297 1205730 1208537 1208615 1210048 1213104 1215017 1215442 1219758 1222577 1222995 1224815 1226234 1227135 1229057 1229469 1231215 1232706 1233609 1234865 1237695 1239004 1241336 1243175 1245255 1245655 1246335 1247810 1251047 1251956 1253206 1253895 1254975 1256442 1257521 1258432 1262535 1262615 1264695 1268147 1268315 1269390 1272255 1273709 1274621 1275495 1276946 1278024 1278813 1280895 1281095 1282341 1287570 1288455 1289894 1291695 1293130 1296204 1296367 1297095 1297173 1298526 1299255 1301762 1302679 1306815 1307159 1308236 1308975 1310232 1311473 1313629 1315454 1316613 1317614 1317948 1319024 1319945 1324423 1325174 1325498 1326419 1329573 1329815 1330736 1330894 1331654 1333050 1336287 1337211 1337366 1339214 1340373 1341374 1341682 1343685 1344918 1348934 1350313 1351391 1352253 1353399 1354628 1355414 1356786 1359872 1361097 1364054 1364133 1367426 1368651 1369454 1370665 1371888 1372773 1376201 1377135 1378094 1379294 1379438 1380335 1384574 1385767 1386992 1387896 1387921 1389000 1389149 1391138 1391152 1392219 1392222 1392225 1392226 1392229 1392231 1394540 1395372 1396699 1397776 1403174 1405092 1405330 1406406 1412880 1415892 1417197 1418273 1422372 1423672 1424747 1429064 1432302 1433172 1434459 1438777 1439652 1439853 1445252 1446327 1450646 1451529 1453882 1456929 1460356 1464675 1466628 1466629 1466830 1467706 1467709 1468809 1468882 1468885 1468886 1469865 1469872 1470952 1471034 1471039 1471049 1471148 1472113 1472131 1473102 1473112 1474180 1474270 1474292 1475289 1475348 1475372 1476353 1476454 1476543 1477433 1479577 1479593 1479664 1479781 1481855 1483020 1483915 1483929 1483982 1484972 1484995 1485061 1486258 1487221 1487255 1488211 1488235 1489329 1490395 1490459 1491539 1492731 1493635 1493699 1494715 1494815 1495809 1495970 1496847 1496939 1497956 1499005 1499049 1499135 1500214 1500288 1501196 1502447 1503417 1504438 1505529 1507642 1507678 1507736 1508816 1508852 1508921 1509932 1510918 1512009 1512159 1513040 1513172 1514158 1515249 1515332 1516318 1517452 1518572 1518634 1519558 1521729 1521772 1522756 1522809 1522954 1523834 1523878 1525012 1527118 1527170 1527210 1528209 1528250 1529232 1529278 1529370 1529431 1531490 1532518 1532529 1532570 1532609 1532669 1533598 1534678 1535769 1536849 1537868 1537918 1537969 1539145 1540026 1540078 1541209 1541249 1544398 1545478 1545489 1545623 1546505 1546569 1547727 1548805 1549743 1549809 1549847 1550925 1551958 1552100 1553049 1555198 1555209 1556219 1556365 1557401 1558449 1558481 1559458 1559529 1560598 1561815 1562695 1562842 1563774 1563849 1567014 1567089 1567117 1567120 1567213 1568092 1568156 1569249 1569321 1570356 1572411 1572561 1573556 1573594 1573598 1573641 1573689 1575649 1575729 1576927 1577914 1578887 1579039 1581129 1582196 1582232 1583311 1583316 1583359 1584286 1584369 1584486 1586554 1587522 1587629 1587679 1588676 1588758 1589682 1589769 1590867 1591916 1591929 1591947 1591998 1592921 1593122 1595169 1595192 1596159 1596249 1596265 1597316 1598432 1598476 1598517 1600569 1600585 1601664 1602752 1602794 1603714 1603809 1603873 1604904 1605956 1606952 1607036 1607113 1607155 1608234 1609112 1609222 1609271 1610289 1612351 1612460 1612510 1613429 1613529 1613588 1614596 1614609 1615676 1615700 1615747 1615790 1616668 1616780 1617849 1617950 1618827 1620064 1622169 1622180 1622268 1623145 1623234 1623302 1624224 1624329 1625420 1626474 1626500 1627621 1628542 1629621 1629740 1629827 1630809 1630859 1630905 1631874 1631988 1631989 1632858 1633067 1634060 1635129 1635177 1635225 1636194 1637174 1637300 1638416 1638464 1640576 1640622 1641491 1641609 1641620 1642674 1642700 1642780 1643858 1644849 1644894 1644944 1645809 1645914 1645940 1645974 1648100 1648177 1649046 1649213 1650345 1651204 1651314 1651329 1651340 1651413 1652492 1653531 1654440 1654569 1654580 1654651 1655519 1655660 1655731 1656771 1656809 1658755 1658874 1659980 1660049 1661034 1661049 1661089 1661128 1663288 1664154 1664302 1664368 1665382 1666310 1667512 1667569 1668622 1668686 1669549 1669689 1670769 1670809 1671705 1671889 1671926 1672943 1674941 1675072 1675128 1675164 1676152 1676183 1677097 1679409 1680562 1681415 1681585 1683572 1683745 1683802 1684792 1684809 1684848 1684882 1685928 1685960 1686985 1687888 1688049 1688065 1688967 1689168 1689200 1691272 1691360 1692203 1692408 1693449 1694359 1694512 1694600 1695627 1695680 1696707 1697596 1699965 1700834 1700992 1701009 1701027 1701078 1703187 1703205 1704249 1704318 1706226 1706427 1706445 1707558 1708552 1709632 1709649 1710545 1710798 1711621 1711827 1712925 1714005 1715067 1715117 1715938 1716112 1717016 1717245 1719369 1719405 1719437 1720432 1720449 1720467 1720517 1722645 1723757 1724569 1724805 1725647 1725832 1725867 1725917 1727803 1728009 1728045 1729072 1729107 1729125 1729157 1730237 1732118 1732347 1732365 1733194 1733409 1733477 1734507 1735351 1735552 1735637 1738827 1739957 1740744 1741821 1742032 1742067 1742085 1743129 1743147 1743197 1743979 1745272 1745325 1746387 1746437 1747215 1747517 1748529 1748565 1749609 1751527 1751752 1751837 1752867 1754761 1754990 1756107 1756157 1756918 1757169 1757187 1757205 1757996 1759397 1761232 1761507 1762311 1762637 1763630 1763649 1763685 1763718 1764464 1764468 1764710 1765543 1766958 1767700 1767950 1768774 1768778 1769067 1769085 1769120 1770129 1770149 1770933 1771190 1772007 1772010 1772360 1773369 1773441 1774165 1774310 1774312 1774315 1774320 1774327 1774330 1774335 1774342 1774345 1774352 1774355 1774469 1775239 1775242 1775510 1776466 1776467 1776520 1776523 1776684 1777394 1777542 1777604 1778805 1779830 1779849 1779925 1780631 1780778 1780847 1780910 1780949 1781929 1782786 1782931 1783089 1783169 1783861 1783862 1784009 1784090 1784150 1784191 1784205 1786015 1786019 1786166 1786365 1787089 1787092 1787332 1787390 1787490 1788511 1788574 1789244 1789247 1789495 1789550 1789569 1790477 1790685 1791398 1791401 1792475 1792632 1792790 1792809 1792831 1792895 1793979 1794630 1794632 1794991 1795005 1795704 1795707 1795865 1796030 1796085 1798939 1798943 1799102 1799220 1799289 1799311 1799382 1800350 1800407 1801095 1801545 1802172 1802464 1802510 1802567 1803416 1804326 1804713 1804788 1805750 1805809 1805870 1805871 1806480 1806482 1806830 1807557 1807730 1807929 1808032 1808806 1808948 1809033 1809049 1809114 1809712 1810070 1811111 1811150 1811277 1811867 1811868 1812289 1812358 1812941 1812943 1813117 1813120 1813329 1814351 1814390 1814433 1814450 1814523 1815095 1815098 1815604 1816174 1816550 1817432 1818753 1818770 1818847 1819409 1819412 1819809 1819930 1820667 1821565 1821746 1822012 1822092 1822643 1823030 1823175 1823903 1824076 1824797 1824800 1824977 1825235 1825340 1825874 1826270 1827134 1827318 1828028 1828032 1828451 1828494 1828583 1829105 1829575 1830558 1830635 1830744 1831263 1831445 1831828 1832337 1832340 1832750 1832817 1833599 1833798 1833851 1833875 1834494 1835568 1835572 1835958 1836058 1838195 1838314 1838805 1838989 1838994 1839230 1840278 1840477 1840960 1840961 1841358 1841411 1841458 1841560 1842037 1842224 1842540 1843721 1844191 1844193 1844380 1844596 1845781 1846348 1846534 1846756 1846790 1846836 1846965 1847423 1847424 1847612 1847870 1848971 1849577 1849580 1849766 1850078 1850103 1850655 1850844 1851288 1852808 1852811 1853001 1853452 1853880 1853887 1854316 1854398 1854536 1855394 1855451 1855506 1856031 1856036 1856039 1856510 1857611 1858391 1858395 1858720 1859263 1859266 1859470 1859714 1859938 1859941 1860338 1860341 1860794 1861022 1861627 1861960 1862492 1862703 1862990 1863186 1864121 1864269 1864644 1864645 1864648 1864860 1865171 1865201 1865226 1865716 1865721 1866253 1866307 1866430 1866434 1867014 1867516 1867868 1867871 1868352 1868441 1868939 1868946 1869170 1869432 1869493 1869681 1870550 1870644 1871095 1871325 1871630 1871681 1871724 1871725 1871727 1871728 1871845 1871848 1872401 1872929 1873248 1873250 1873867 1873890 1873891 1873893 1874872 1874893 1874949 1874964 1874974 1875093 1875095 1875399 1875401 1875910 1877083 1877553 1877556 1877791 1877795 1878133 1878163 1878215 1878342 1878870 1879284 1879297 1879423 1879426 1879707 1880378 1880380 1881373 1881405 1881429 1881445 1881593 1881596 1881863 1882453 1882542 1882683 1883185 1883512 1883623 1884017 1884257 1884261 1884549 1884592 1884647 1884671 1884685 1884844 1884847 1885094 1885627 1885693 1885785 1885931 1885933 1886416 1886807 1886831 1886866 1887244 1887249 1887493 1887854 1887927 1888096 1888099 1888912 1889029 1889031 1889393 1889399 1889651 1890071 1890087 1890263 1890269 1890725 1890728 1891094 1891151 1891350 1891353 1891358 1891541 1891545 1891549 1891550 1892154 1892209 1892276 1893290 1893328 1893357 1893361 1893519 1893527 1893528 1893531 1893693 1893696 1894265 1894334 1894392 1894764 1894767 1894768 1895035 1895039 1895394 1895452 1895525 1895697 1895700 1895701 1895705 1895711 1896474 1897613 1897632 1897648 1897686 1897689 1897874 1897875 1897883 1897884 1897978 1897981 1897982 1897985 1897988 1897990 1897994 1897995 1897998 1898272 1898654 1898694 1898712 1898730 1899347 1899350 1899850 1900051 1900052 1900056 1900062 1900065 1900068 1900073 1900075 1900082 1900084 1900091 1900094 1900095 1900104 1900107 1900108 1900111 1900113 1900117 1900122 1900126 1900129 1900130 1900134 1900814 1900872 1900890 1901503 1901874 1901949 1901951 1902014 1902017 1902902 1903018 1903028 1903099 1903653 1903654 1903658 1903980 1904054 1904131 1904725 1905114 1905183 1905184 1905186 1905262 1905264 1906261 1906349 1906868 1906874 1906877 1906884 1907220 1907294 1907373 1908354 1908513 1908517 1908913 1908916 1908922 1908927 1908931 1908936 1908944 1908945 1908950 1908954 1908960 1908963 1908968 1908972 1908978 1908983 1908987 1908993 1908998 1909006 1909012 1909019 1909020 1909434 1910534 1910680 1910684 1911058 1911062 1911067 1911539 1911614 1911695 1912129 1912132 1912674 1912849 1913932 1913935 1913938 1914279 1914283 1914777 1914834 1914854 1914937 1915857 1917181 1917514 1918177 1919154 1919346 1919348 1919351 1919662 1919666 1920254 1920437 1921255 1921420 1921801 1921805 1921811 1921816 1922335 1922394 1922601 1922604 1923474 1923943 1923946 1923953 1923956 1924662 1924769 1925012 1925013 1925019 1925022 1925654 1925853 1925855 1927153 1927156 1927157 1927160 1927163 1927165 1927168 1928218 1928220 1928223 1928228 1928811 1928894 1928986 1930067 1930182 1930183 1930362 1930366 1930371 1930375 1931034 1931054 1932051 1932350 1932353 1932504 1932507 1932512 1932515 1933311 1934294 1934640 1934642 1934645 1934650 1934653 1934655 1934659 1936456 1936552 1936679 1936682 1936768 1936771 1936776 1936780 1936783 1936785 1936791 1936793 1936796 1938527 1938849 1938851 1938854 1938859 1938862 1938864 1938871 1938873 1938878 1938881 1938888 1938891 1938896 1938904 1938905 1938913 1938918 1938919 1939674 1939696 1940776 1942847 1944016 1944120 1945074 1947365 1948243 1948336 1949529 1951554 1951576 1952562 1953852 1954816 1958034 1959040 1959136 1959261 1960347 1961296 1962277 1962510 1964677 1965616 1966674 1968751 1969002 1969936 1970086 1971016 1971726 1971733 1971744 1971750 1971761 1971766 1971771 1971780 1971785 1971790 1971792 1971799 1971803 1971808 1972074 1972251 1972491 1972495 1972501 1972506 1972509 1972520 1972525 1972534 1972539 1972550 1972556 1972563 1972574 1972579 1972590 1972596 1973681 1973688 1973693 1973696 1973702 1973706 1973711 1973718 1973841 1973846 1973853 1973863 1973869 1973874 1973973 1973978 1974256 1974418 1974626 1974632 1974638 1974643 1975230 1975314 1976772 1976781 1977223 1977226 1977228 1977496 1978050 1978061 1978071 1978078 1978088 1978099 1978106 1978117 1978128 1978139 1978150 1978916 1978922 1979395 1979546 1979826 1979827 1980736 1981063 1981070 1981563 1981568 1981794 1981816 1981990 1982655 1982658 1983216 1983862 1984157 1984160 1984165 1984824 1984829 1985056 1985249 1985255 1985368 1985370 1986018 1986998 1987419 1987426 1987514 1987520 1988276 1988296 1988571 1988573 1988579 1988585 1988591 1989159 1989164 1989167 1989169 1989173 1989174 1989177 1989180 1989254 1989592 1989598 1989601 1989606 1989609 1990269 1990274 1990277 1990283 1990289 1990297 1990298 1990692 1990697 1990700 1990707 1990710 1990711 1990718 1990719 1990722 1991411 1991536 1992464 1992467 1992472 1992476 1992477 1992480 1992481 1992485 1992487 1992490 1996916 1996936 2001236 2002336 2004476 2005576 2007716 2010956 2010976 2014196 2016376 2018516 2019616 2020676 2021756 2022856 2024996 2026096 2028256 2029336 2030396 2031496 2032556 2032576 2033656 2036503 2036507 2036510 2036511 2036514 2036517 2036523 2036527 2036528 2036531 2036539 2036543 2036548 2036554 2036557 2036558 2036566 2036567 2036570 2036896 2038491 2038495 2038501 2038506 2038510 2038515 2038518 2038524 2038532 2038533 2038539 2038542 2038550 2038556 2038561 2038569 2038575 2038581 2038591 2038595 2038600 2038605 2038612 2038615 2038620 2038624 2038630 2038637 2038639 2038642 2038648 2038650 2038653 2038656 2038661 2038736 2038740 2038741 2038745 2038750 2039036 2040136 2040618 2040623 2040632 2040633 2040640 2040646 2040913 2040916 2040919 2040922 2040925 2040930 2041217 2041674 2041680 2041683 2041689 2041692 2042276 2043094 2043098 2043103 2043822 2043829 2044187 2044192 2044457 2045836 2045843 2045846 2045851 2045853 2045856 2045864 2045870 2045871 2045874 2045878 2045884 2045887 2045892 2045895 2045902 2045909 2045975 2046357 2046360 2046617 2047979 2047981 2047988 2047990 2048076 2048081 2048091 2048095 2048101 2048106 2048114 2048120 2048125 2048128 2048525 2048531 2049857 2050117 2050122 2050126 2050134 2050695 2050696 2050703 2050704 2050916 2050937 2051185 2051191 2051788 2051789 2053099 2053325 2053329 2053334 2053337 2053953 2053956 2053958 2054154 2056121 2056122 2056125 2056547 2056551 2056556 2056559 2057420 2058287 2058291 2058294 2058297 2058299 2058302 2058303 2058306 2058310 2058312 2058315 2058318 2058319 2058322 2058327 2058330 2058331 2058334 2058337 2058502 2058684 2058692 2058697 2058699 2058705 2059421 2059422 2059430 2059431 2059438 2059439 2059731 2059735 2059740 2059744 2059751 2059757 2059763 2061603 2061605 2061608 2061611 2061614 2061617 2061622 2061745 2061746 2061752 2061757 2061869 2061873 2061879 2061882 2061885 2062794 2063787 2063790 2063874 2063919 2063922 2063923 2063925 2063929 2063930 2063933 2063936 2063938 2063941 2063944 2063948 2063951 2063953 2063956 2063959 2063962 2063969 2063972 2063973 2063976 2063979 2063982 2063985 2063988 2063995 2063997 2064003 2064006 2064011 2064014 2064019 2064024 2064874 2064877 2064881 2064883 2064886 2064887 2064890 2064891 2064894 2064897 2064900 2064903 2066034 2067067 2067068 2067071 2067074 2067077 2067080 2067081 2067084 2067085 2067087 2069248 2069251 2069256 2069259 2069260 2069273 2069274 2070342 2070343 2070344 2070346 2070347 2070350 2070353 ];
L = superpixels(X,10468);
BW = lazysnapping(X,L,foregroundInd,backgroundInd);

% Create masked image.
maskedImage = X;
maskedImage(~BW) = 0;
end

