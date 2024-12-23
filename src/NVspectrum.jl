module NVspectrum
const wavelength = [460.7, 460.96, 461.22, 461.49, 461.75, 462.02, 462.28, 462.54, 462.81, 463.07, 
463.34, 463.6, 463.86, 464.13, 464.39, 464.66, 464.92, 465.18, 465.45, 465.71, 
465.98, 466.24, 466.5, 466.77, 467.03, 467.3, 467.56, 467.82, 468.09, 468.35, 
468.62, 468.88, 469.14, 469.41, 469.67, 469.94, 470.2, 470.46, 470.73, 470.99, 
471.26, 471.52, 471.78, 472.05, 472.31, 472.58, 472.84, 473.1, 473.37, 473.63, 
473.9, 474.16, 474.42, 474.69, 474.95, 475.22, 475.48, 475.74, 476.01, 476.27, 
476.54, 476.8, 477.06, 477.33, 477.59, 477.86, 478.12, 478.38, 478.65, 478.91, 
479.18, 479.44, 479.7, 479.97, 480.23, 480.5, 480.76, 481.02, 481.29, 481.55, 
481.81, 482.08, 482.34, 482.61, 482.87, 483.13, 483.4, 483.66, 483.93, 484.19, 
484.45, 484.72, 484.98, 485.25, 485.51, 485.77, 486.04, 486.3, 486.57, 486.83, 
487.09, 487.36, 487.62, 487.89, 488.15, 488.41, 488.68, 488.94, 489.2, 489.47, 
489.73, 490.0, 490.26, 490.52, 490.79, 491.05, 491.32, 491.58, 491.84, 492.11, 
492.37, 492.64, 492.9, 493.16, 493.43, 493.69, 493.96, 494.22, 494.48, 494.75, 
495.01, 495.27, 495.54, 495.8, 496.07, 496.33, 496.59, 496.86, 497.12, 497.39, 
497.65, 497.91, 498.18, 498.44, 498.7, 498.97, 499.23, 499.5, 499.76, 500.02, 
500.29, 500.55, 500.82, 501.08, 501.34, 501.61, 501.87, 502.13, 502.4, 502.66, 
502.93, 503.19, 503.45, 503.72, 503.98, 504.25, 504.51, 504.77, 505.04, 505.3, 
505.56, 505.83, 506.09, 506.36, 506.62, 506.88, 507.15, 507.41, 507.67, 507.94, 
508.2, 508.47, 508.73, 508.99, 509.26, 509.52, 509.79, 510.05, 510.31, 510.58, 
510.84, 511.1, 511.37, 511.63, 511.9, 512.16, 512.42, 512.69, 512.95, 513.21, 
513.48, 513.74, 514.01, 514.27, 514.53, 514.8, 515.06, 515.32, 515.59, 515.85, 
516.12, 516.38, 516.64, 516.91, 517.17, 517.43, 517.7, 517.96, 518.23, 518.49, 
518.75, 519.02, 519.28, 519.54, 519.81, 520.07, 520.34, 520.6, 520.86, 521.13, 
521.39, 521.65, 521.92, 522.18, 522.45, 522.71, 522.97, 523.24, 523.5, 523.76, 
524.03, 524.29, 524.56, 524.82, 525.08, 525.35, 525.61, 525.87, 526.14, 526.4, 
526.67, 526.93, 527.19, 527.46, 527.72, 527.98, 528.25, 528.51, 528.77, 529.04, 
529.3, 529.57, 529.83, 530.09, 530.36, 530.62, 530.88, 531.15, 531.41, 531.68, 
531.94, 532.2, 532.47, 532.73, 532.99, 533.26, 533.52, 533.78, 534.05, 534.31, 
534.58, 534.84, 535.1, 535.37, 535.63, 535.89, 536.16, 536.42, 536.68, 536.95, 
537.21, 537.48, 537.74, 538.0, 538.27, 538.53, 538.79, 539.06, 539.32, 539.58, 
539.85, 540.11, 540.38, 540.64, 540.9, 541.17, 541.43, 541.69, 541.96, 542.22, 
542.48, 542.75, 543.01, 543.27, 543.54, 543.8, 544.07, 544.33, 544.59, 544.86, 
545.12, 545.38, 545.65, 545.91, 546.17, 546.44, 546.7, 546.96, 547.23, 547.49, 
547.76, 548.02, 548.28, 548.55, 548.81, 549.07, 549.34, 549.6, 549.86, 550.13, 
550.39, 550.65, 550.92, 551.18, 551.45, 551.71, 551.97, 552.24, 552.5, 552.76, 
553.03, 553.29, 553.55, 553.82, 554.08, 554.34, 554.61, 554.87, 555.13, 555.4, 
555.66, 555.93, 556.19, 556.45, 556.72, 556.98, 557.24, 557.51, 557.77, 558.03, 
558.3, 558.56, 558.82, 559.09, 559.35, 559.61, 559.88, 560.14, 560.4, 560.67, 
560.93, 561.2, 561.46, 561.72, 561.99, 562.25, 562.51, 562.78, 563.04, 563.3, 
563.57, 563.83, 564.09, 564.36, 564.62, 564.88, 565.15, 565.41, 565.67, 565.94, 
566.2, 566.46, 566.73, 566.99, 567.25, 567.52, 567.78, 568.05, 568.31, 568.57, 
568.84, 569.1, 569.36, 569.63, 569.89, 570.15, 570.42, 570.68, 570.94, 571.21, 
571.47, 571.73, 572.0, 572.26, 572.52, 572.79, 573.05, 573.31, 573.58, 573.84, 
574.1, 574.37, 574.63, 574.89, 575.16, 575.42, 575.68, 575.95, 576.21, 576.47, 
576.74, 577.0, 577.26, 577.53, 577.79, 578.05, 578.32, 578.58, 578.84, 579.11, 
579.37, 579.63, 579.9, 580.16, 580.42, 580.69, 580.95, 581.21, 581.48, 581.74, 
582.0, 582.27, 582.53, 582.79, 583.06, 583.32, 583.58, 583.85, 584.11, 584.37, 
584.64, 584.9, 585.16, 585.43, 585.69, 585.95, 586.22, 586.48, 586.74, 587.01, 
587.27, 587.53, 587.8, 588.06, 588.32, 588.59, 588.85, 589.11, 589.38, 589.64, 
589.9, 590.17, 590.43, 590.69, 590.96, 591.22, 591.48, 591.75, 592.01, 592.27, 
592.54, 592.8, 593.06, 593.33, 593.59, 593.85, 594.12, 594.38, 594.64, 594.91, 
595.17, 595.43, 595.7, 595.96, 596.22, 596.49, 596.75, 597.01, 597.28, 597.54, 
597.8, 598.06, 598.33, 598.59, 598.85, 599.12, 599.38, 599.64, 599.91, 600.17, 
600.43, 600.7, 600.96, 601.22, 601.49, 601.75, 602.01, 602.28, 602.54, 602.8, 
603.07, 603.33, 603.59, 603.86, 604.12, 604.38, 604.64, 604.91, 605.17, 605.43, 
605.7, 605.96, 606.22, 606.49, 606.75, 607.01, 607.28, 607.54, 607.8, 608.07, 
608.33, 608.59, 608.86, 609.12, 609.38, 609.64, 609.91, 610.17, 610.43, 610.7, 
610.96, 611.22, 611.49, 611.75, 612.01, 612.28, 612.54, 612.8, 613.07, 613.33, 
613.59, 613.85, 614.12, 614.38, 614.64, 614.91, 615.17, 615.43, 615.7, 615.96, 
616.22, 616.49, 616.75, 617.01, 617.27, 617.54, 617.8, 618.06, 618.33, 618.59, 
618.85, 619.12, 619.38, 619.64, 619.91, 620.17, 620.43, 620.69, 620.96, 621.22, 
621.48, 621.75, 622.01, 622.27, 622.54, 622.8, 623.06, 623.33, 623.59, 623.85, 
624.11, 624.38, 624.64, 624.9, 625.17, 625.43, 625.69, 625.96, 626.22, 626.48, 
626.74, 627.01, 627.27, 627.53, 627.8, 628.06, 628.32, 628.59, 628.85, 629.11, 
629.37, 629.64, 629.9, 630.16, 630.43, 630.69, 630.95, 631.22, 631.48, 631.74, 
632.0, 632.27, 632.53, 632.79, 633.06, 633.32, 633.58, 633.84, 634.11, 634.37, 
634.63, 634.9, 635.16, 635.42, 635.69, 635.95, 636.21, 636.47, 636.74, 637.0, 
637.26, 637.53, 637.79, 638.05, 638.31, 638.58, 638.84, 639.1, 639.37, 639.63, 
639.89, 640.15, 640.42, 640.68, 640.94, 641.21, 641.47, 641.73, 642.0, 642.26, 
642.52, 642.78, 643.05, 643.31, 643.57, 643.84, 644.1, 644.36, 644.62, 644.89, 
645.15, 645.41, 645.68, 645.94, 646.2, 646.46, 646.73, 646.99, 647.25, 647.52, 
647.78, 648.04, 648.3, 648.57, 648.83, 649.09, 649.36, 649.62, 649.88, 650.14, 
650.41, 650.67, 650.93, 651.19, 651.46, 651.72, 651.98, 652.25, 652.51, 652.77, 
653.03, 653.3, 653.56, 653.82, 654.09, 654.35, 654.61, 654.87, 655.14, 655.4, 
655.66, 655.92, 656.19, 656.45, 656.71, 656.98, 657.24, 657.5, 657.76, 658.03, 
658.29, 658.55, 658.82, 659.08, 659.34, 659.6, 659.87, 660.13, 660.39, 660.65, 
660.92, 661.18, 661.44, 661.71, 661.97, 662.23, 662.49, 662.76, 663.02, 663.28, 
663.54, 663.81, 664.07, 664.33, 664.6, 664.86, 665.12, 665.38, 665.65, 665.91, 
666.17, 666.43, 666.7, 666.96, 667.22, 667.48, 667.75, 668.01, 668.27, 668.54, 
668.8, 669.06, 669.32, 669.59, 669.85, 670.11, 670.37, 670.64, 670.9, 671.16, 
671.42, 671.69, 671.95, 672.21, 672.47, 672.74, 673.0, 673.26, 673.53, 673.79, 
674.05, 674.31, 674.58, 674.84, 675.1, 675.36, 675.63, 675.89, 676.15, 676.41, 
676.68, 676.94, 677.2, 677.46, 677.73, 677.99, 678.25, 678.51, 678.78, 679.04, 
679.3, 679.57, 679.83, 680.09, 680.35, 680.62, 680.88, 681.14, 681.4, 681.67, 
681.93, 682.19, 682.45, 682.72, 682.98, 683.24, 683.5, 683.77, 684.03, 684.29, 
684.55, 684.82, 685.08, 685.34, 685.6, 685.87, 686.13, 686.39, 686.65, 686.92, 
687.18, 687.44, 687.7, 687.97, 688.23, 688.49, 688.75, 689.02, 689.28, 689.54, 
689.8, 690.07, 690.33, 690.59, 690.85, 691.12, 691.38, 691.64, 691.9, 692.17, 
692.43, 692.69, 692.95, 693.22, 693.48, 693.74, 694.0, 694.26, 694.53, 694.79, 
695.05, 695.31, 695.58, 695.84, 696.1, 696.36, 696.63, 696.89, 697.15, 697.41, 
697.68, 697.94, 698.2, 698.46, 698.73, 698.99, 699.25, 699.51, 699.78, 700.04, 
700.3, 700.56, 700.83, 701.09, 701.35, 701.61, 701.87, 702.14, 702.4, 702.66, 
702.92, 703.19, 703.45, 703.71, 703.97, 704.24, 704.5, 704.76, 705.02, 705.29, 
705.55, 705.81, 706.07, 706.33, 706.6, 706.86, 707.12, 707.38, 707.65, 707.91, 
708.17, 708.43, 708.7, 708.96, 709.22, 709.48, 709.74, 710.01, 710.27, 710.53, 
710.79, 711.06, 711.32, 711.58, 711.84, 712.1, 712.37, 712.63, 712.89, 713.15, 
713.42, 713.68, 713.94, 714.2, 714.47, 714.73, 714.99, 715.25, 715.51, 715.78, 
716.04, 716.3, 716.56, 716.83, 717.09, 717.35, 717.61, 717.87, 718.14, 718.4, 
718.66, 718.92, 719.19, 719.45, 719.71, 719.97, 720.23, 720.5, 720.76, 721.02, 
721.28, 721.54, 721.81, 722.07, 722.33, 722.59, 722.86, 723.12, 723.38, 723.64, 
723.9, 724.17, 724.43, 724.69, 724.95, 725.22, 725.48, 725.74, 726.0, 726.26, 
726.53, 726.79, 727.05, 727.31, 727.57, 727.84, 728.1, 728.36, 728.62, 728.88, 
729.15, 729.41, 729.67, 729.93, 730.2, 730.46, 730.72, 730.98, 731.24, 731.51, 
731.77, 732.03, 732.29, 732.55, 732.82, 733.08, 733.34, 733.6, 733.86, 734.13, 
734.39, 734.65, 734.91, 735.17, 735.44, 735.7, 735.96, 736.22, 736.48, 736.75, 
737.01, 737.27, 737.53, 737.79, 738.06, 738.32, 738.58, 738.84, 739.1, 739.37, 
739.63, 739.89, 740.15, 740.41, 740.68, 740.94, 741.2, 741.46, 741.72, 741.99, 
742.25, 742.51, 742.77, 743.03, 743.3, 743.56, 743.82, 744.08, 744.34, 744.61, 
744.87, 745.13, 745.39, 745.65, 745.92, 746.18, 746.44, 746.7, 746.96, 747.23, 
747.49, 747.75, 748.01, 748.27, 748.54, 748.8, 749.06, 749.32, 749.58, 749.84, 
750.11, 750.37, 750.63, 750.89, 751.15, 751.42, 751.68, 751.94, 752.2, 752.46, 
752.73, 752.99, 753.25, 753.51, 753.77, 754.03, 754.3, 754.56, 754.82, 755.08, 
755.34, 755.61, 755.87, 756.13, 756.39, 756.65, 756.92, 757.18, 757.44, 757.7, 
757.96, 758.22, 758.49, 758.75, 759.01, 759.27, 759.53, 759.8, 760.06, 760.32, 
760.58, 760.84, 761.1, 761.37, 761.63, 761.89, 762.15, 762.41, 762.67, 762.94, 
763.2, 763.46, 763.72, 763.98, 764.25, 764.51, 764.77, 765.03, 765.29, 765.55, 
765.82, 766.08, 766.34, 766.6, 766.86, 767.12, 767.39, 767.65, 767.91, 768.17, 
768.43, 768.69, 768.96, 769.22, 769.48, 769.74, 770.0, 770.26, 770.53, 770.79, 
771.05, 771.31, 771.57, 771.83, 772.1, 772.36, 772.62, 772.88, 773.14, 773.4, 
773.67, 773.93, 774.19, 774.45, 774.71, 774.97, 775.24, 775.5, 775.76, 776.02, 
776.28, 776.54, 776.81, 777.07, 777.33, 777.59, 777.85, 778.11, 778.38, 778.64, 
778.9, 779.16, 779.42, 779.68, 779.95, 780.21, 780.47, 780.73, 780.99, 781.25, 
781.51, 781.78, 782.04, 782.3, 782.56, 782.82, 783.08, 783.35, 783.61, 783.87, 
784.13, 784.39, 784.65, 784.91, 785.18, 785.44, 785.7, 785.96, 786.22, 786.48, 
786.75, 787.01, 787.27, 787.53, 787.79, 788.05, 788.31, 788.58, 788.84, 789.1, 
789.36, 789.62, 789.88, 790.14, 790.41, 790.67, 790.93, 791.19, 791.45, 791.71, 
791.98, 792.24, 792.5, 792.76, 793.02, 793.28, 793.54, 793.81, 794.07, 794.33, 
794.59, 794.85, 795.11, 795.37, 795.64, 795.9, 796.16, 796.42, 796.68, 796.94, 
797.2, 797.46, 797.73, 797.99, 798.25, 798.51, 798.77, 799.03, 799.29, 799.56, 
799.82, 800.08, 800.34, 800.6, 800.86, 801.12, 801.39, 801.65, 801.91, 802.17, 
802.43, 802.69, 802.95, 803.21, 803.48, 803.74, 804.0, 804.26, 804.52, 804.78, 
805.04, 805.31, 805.57, 805.83, 806.09, 806.35, 806.61, 806.87, 807.13, 807.4, 
807.66, 807.92, 808.18, 808.44, 808.7, 808.96, 809.22, 809.49, 809.75, 810.01, 
810.27, 810.53, 810.79, 811.05, 811.31, 811.58, 811.84, 812.1, 812.36, 812.62, 
] * 1e-3

const spectrum = [-206.0, -70.0, -67.0, 40.0, -28.0, 55.0, 29.0, 80.0, -61.0, 35.0, 
21.0, -62.0, -77.0, -18.0, 54.0, -90.0, -47.0, -26.0, -53.0, -4.0, 
41.0, 28.0, 15.0, -10.0, -22.0, 47.0, -67.0, 9.0, 13.0, -12.0, 
-8.0, 6.0, -56.0, 44.0, 45.0, 56.0, -83.0, 96.0, -65.0, 44.0, 
-7.0, 35.0, -23.0, 14.0, -14.0, 36.0, -51.0, -42.0, -31.0, 10.0, 
10.0, 25.0, -58.0, -74.0, 76.0, -118.0, 9.0, -56.0, 1.0, -87.0, 
-24.0, 16.0, -97.0, 47.0, -62.0, 3.0, -40.0, 19.0, 67.0, -162.0, 
61.0, 23.0, 94.0, -16.0, -19.0, -19.0, 40.0, 22.0, -77.0, -5.0, 
2.0, 62.0, 67.0, -41.0, -24.0, -88.0, 29.0, 4.0, 42.0, 39.0, 
-46.0, 107.0, 28.0, 22.0, 87.0, -61.0, 97.0, -3.0, 43.0, 0.0, 
-57.0, -7.0, -80.0, 25.0, 2.0, -82.0, -8.0, 106.0, 90.0, 8.0, 
-41.0, -15.0, -6.0, -29.0, -9.0, -36.0, 36.0, 20.0, -33.0, 51.0, 
-119.0, -1.0, -19.0, 24.0, -167.0, -45.0, -112.0, -30.0, -92.0, -120.0, 
60.0, 63.0, -68.0, 5.0, -61.0, -23.0, -12.0, 11.0, 15.0, 34.0, 
-38.0, 19.0, -46.0, 13.0, -78.0, 28.0, 35.0, -37.0, -41.0, -6.0, 
13.0, 38.0, 1.0, 17.0, -14.0, 8.0, 15.0, -62.0, 10.0, 99.0, 
108.0, -31.0, 5.0, 30.0, 92.0, -111.0, 120.0, -24.0, -10.0, 36.0, 
-38.0, 17.0, -26.0, 74.0, -56.0, -95.0, 6.0, 6.0, 26.0, 23.0, 
-3.0, 36.0, 58.0, -75.0, 77.0, -23.0, -129.0, -43.0, 10.0, -24.0, 
78.0, 31.0, -26.0, -38.0, 21.0, -24.0, 178.0, 248.0, 340.0, 482.0, 
475.0, 442.0, 590.0, 492.0, 518.0, 496.0, 409.0, 274.0, 132.0, 41.0, 
47.0, 33.0, 94.0, 57.0, 43.0, 126.0, -35.0, 84.0, 102.0, -15.0, 
-86.0, -88.0, -14.0, 23.0, -30.0, 66.0, 33.0, -17.0, 110.0, -27.0, 
34.0, -65.0, 20.0, 18.0, -9.0, -12.0, 6.0, 32.0, -19.0, 54.0, 
42.0, -35.0, -30.0, -43.0, 103.0, 143.0, -51.0, -39.0, -7.0, -80.0, 
-48.0, -16.0, 68.0, 34.0, -20.0, 43.0, 78.0, -14.0, -48.0, 124.0, 
28.0, 70.0, -14.0, 42.0, -2.0, -87.0, -26.0, -77.0, 55.0, -73.0, 
9.0, -4.0, -50.0, 76.0, -59.0, 74.0, -19.0, 103.0, 13.0, -8.0, 
58.0, 45.0, 82.0, -62.0, 9.0, -12.0, 45.0, 45.0, 9.0, 138.0, 
-29.0, 41.0, -77.0, 31.0, 4.0, 87.0, 84.0, 61.0, -36.0, 33.0, 
151.0, 124.0, 134.0, 16.0, 149.0, 165.0, 154.0, 173.0, 132.0, 250.0, 
156.0, 264.0, 164.0, 223.0, 324.0, 179.0, 210.0, 211.0, 223.0, 282.0, 
228.0, 168.0, 198.0, 159.0, 58.0, 146.0, 89.0, 136.0, 69.0, 109.0, 
79.0, 19.0, 195.0, 46.0, 162.0, 80.0, 128.0, 62.0, 95.0, 189.0, 
40.0, 50.0, -70.0, 29.0, 67.0, 38.0, -8.0, 27.0, 17.0, 7.0, 
50.0, 27.0, 121.0, 22.0, -25.0, -1.0, 44.0, 103.0, 8.0, 67.0, 
-75.0, 24.0, -22.0, 13.0, 79.0, 40.0, 34.0, 16.0, 31.0, 35.0, 
25.0, 16.0, -2.0, 26.0, -26.0, 48.0, -21.0, 73.0, -7.0, -13.0, 
12.0, -13.0, 16.0, -36.0, -17.0, 45.0, 154.0, -34.0, 142.0, 48.0, 
-3.0, 53.0, 37.0, 24.0, 27.0, 63.0, 33.0, 59.0, -19.0, 13.0, 
-18.0, 49.0, 25.0, 52.0, -104.0, 86.0, -37.0, 18.0, 80.0, -43.0, 
45.0, 25.0, 62.0, 29.0, 89.0, 111.0, 28.0, 111.0, -73.0, -79.0, 
92.0, 100.0, -4.0, -78.0, 104.0, 62.0, 24.0, 51.0, -14.0, -15.0, 
68.0, -52.0, -86.0, 66.0, -5.0, 103.0, -60.0, 63.0, 55.0, 49.0, 
23.0, -5.0, 41.0, 19.0, -87.0, 23.0, 26.0, 71.0, 18.0, 60.0, 
61.0, 86.0, 8.0, 80.0, 69.0, 97.0, 71.0, -32.0, 13.0, 88.0, 
154.0, 140.0, 106.0, 87.0, 1.0, 75.0, 109.0, 114.0, 73.0, -13.0, 
33.0, 200.0, 35.0, 85.0, 109.0, 17.0, 57.0, 132.0, 41.0, 133.0, 
132.0, 57.0, 38.0, 4.0, 13.0, 33.0, 32.0, -5.0, 29.0, 162.0, 
105.0, 140.0, 88.0, 84.0, 19.0, 151.0, 26.0, 10.0, 128.0, 62.0, 
7.0, 85.0, 76.0, 26.0, 64.0, 68.0, 90.0, 82.0, 90.0, 9.0, 
48.0, 199.0, 213.0, 143.0, 299.0, 257.0, 403.0, 369.0, 589.0, 716.0, 
1061.0, 1276.0, 1629.0, 2067.0, 2534.0, 3028.0, 3328.0, 3510.0, 3836.0, 4021.0, 
4111.0, 4259.0, 4133.0, 4188.0, 4126.0, 4111.0, 4023.0, 4031.0, 4196.0, 4025.0, 
4034.0, 4229.0, 4050.0, 4068.0, 4137.0, 4140.0, 4202.0, 4213.0, 4064.0, 4332.0, 
4279.0, 4290.0, 4225.0, 4119.0, 4280.0, 4211.0, 4244.0, 4410.0, 4563.0, 4288.0, 
4504.0, 4525.0, 4631.0, 4904.0, 4811.0, 4998.0, 4971.0, 5189.0, 5124.0, 5213.0, 
5127.0, 5367.0, 5158.0, 5363.0, 5079.0, 5436.0, 5349.0, 5452.0, 5459.0, 5531.0, 
5536.0, 5669.0, 5674.0, 5651.0, 5808.0, 6005.0, 5803.0, 5762.0, 6129.0, 5987.0, 
6282.0, 6111.0, 6199.0, 6435.0, 6329.0, 6330.0, 6745.0, 6424.0, 6640.0, 6613.0, 
6547.0, 6622.0, 6771.0, 6580.0, 6653.0, 6904.0, 6933.0, 6726.0, 6784.0, 7020.0, 
6978.0, 6897.0, 6955.0, 6882.0, 7104.0, 7009.0, 6970.0, 7168.0, 7313.0, 7088.0, 
7310.0, 7275.0, 6958.0, 7227.0, 7282.0, 7192.0, 7408.0, 7389.0, 7431.0, 7318.0, 
7357.0, 7699.0, 7495.0, 7713.0, 7723.0, 7805.0, 7852.0, 8242.0, 8094.0, 7972.0, 
8315.0, 7841.0, 8282.0, 8565.0, 8306.0, 8478.0, 8350.0, 8601.0, 8848.0, 8725.0, 
8900.0, 8992.0, 8872.0, 9409.0, 9525.0, 9782.0, 9705.0, 10113.0, 10836.0, 10828.0, 
11637.0, 12059.0, 12969.0, 13474.0, 14303.0, 14893.0, 14853.0, 15138.0, 14668.0, 13872.0, 
13964.0, 13336.0, 12670.0, 12962.0, 12398.0, 12306.0, 12498.0, 12344.0, 11990.0, 12316.0, 
11828.0, 11991.0, 12519.0, 12171.0, 12305.0, 12318.0, 12327.0, 12414.0, 12358.0, 12530.0, 
12699.0, 12305.0, 12753.0, 12655.0, 12948.0, 12870.0, 12924.0, 13056.0, 13237.0, 13201.0, 
13213.0, 13125.0, 13453.0, 13337.0, 13668.0, 13865.0, 14212.0, 13863.0, 14310.0, 14180.0, 
14132.0, 14286.0, 14724.0, 14650.0, 14952.0, 15151.0, 15432.0, 15229.0, 15754.0, 15722.0, 
15645.0, 16072.0, 16058.0, 16386.0, 16598.0, 16485.0, 16870.0, 16862.0, 17152.0, 17479.0, 
17853.0, 17913.0, 18048.0, 18036.0, 18416.0, 18517.0, 19104.0, 19076.0, 19245.0, 19207.0, 
19467.0, 19714.0, 19779.0, 20019.0, 20066.0, 20376.0, 20858.0, 20677.0, 20734.0, 20526.0, 
20523.0, 20792.0, 20568.0, 20462.0, 20515.0, 20980.0, 20196.0, 21013.0, 20903.0, 20324.0, 
20302.0, 20370.0, 20355.0, 20274.0, 20313.0, 20196.0, 19783.0, 20058.0, 20398.0, 20143.0, 
20395.0, 20331.0, 19941.0, 19961.0, 20223.0, 19825.0, 20041.0, 19609.0, 19643.0, 19733.0, 
19641.0, 20346.0, 19995.0, 20234.0, 20425.0, 20310.0, 20273.0, 20473.0, 20269.0, 19943.0, 
20304.0, 20285.0, 20366.0, 20047.0, 20835.0, 20708.0, 20679.0, 21146.0, 21237.0, 21278.0, 
21551.0, 21426.0, 21006.0, 21463.0, 21320.0, 21132.0, 21213.0, 21603.0, 21650.0, 21958.0, 
22109.0, 22102.0, 22638.0, 22332.0, 22240.0, 22508.0, 22434.0, 22603.0, 22410.0, 22128.0, 
22201.0, 22884.0, 23168.0, 23149.0, 23317.0, 23060.0, 23345.0, 23195.0, 23515.0, 23236.0, 
23207.0, 23085.0, 22781.0, 22775.0, 23561.0, 23647.0, 23379.0, 23665.0, 23819.0, 23631.0, 
23778.0, 23551.0, 23807.0, 23502.0, 23288.0, 23364.0, 23351.0, 23200.0, 22715.0, 23351.0, 
23212.0, 23425.0, 23319.0, 23319.0, 23620.0, 23134.0, 23063.0, 23066.0, 22865.0, 22611.0, 
22691.0, 22566.0, 22428.0, 22669.0, 22322.0, 22761.0, 22514.0, 22491.0, 22220.0, 22042.0, 
21557.0, 21695.0, 21809.0, 22205.0, 21830.0, 21395.0, 21762.0, 21508.0, 21834.0, 22167.0, 
22366.0, 21870.0, 21777.0, 21838.0, 21791.0, 21296.0, 21531.0, 21367.0, 21448.0, 21084.0, 
21197.0, 21272.0, 21296.0, 21897.0, 21871.0, 21495.0, 21665.0, 21702.0, 21633.0, 20979.0, 
21329.0, 20705.0, 20992.0, 20714.0, 20941.0, 21421.0, 21004.0, 21525.0, 21337.0, 21182.0, 
21250.0, 21510.0, 21365.0, 21019.0, 21508.0, 20800.0, 20616.0, 20650.0, 20731.0, 20362.0, 
20875.0, 20798.0, 20625.0, 20983.0, 21146.0, 20700.0, 20849.0, 20742.0, 20851.0, 20945.0, 
20504.0, 20347.0, 20305.0, 20112.0, 20123.0, 19697.0, 19814.0, 19934.0, 20148.0, 20139.0, 
20171.0, 20042.0, 19998.0, 19909.0, 20018.0, 19825.0, 19487.0, 19031.0, 18951.0, 18900.0, 
19109.0, 18765.0, 18985.0, 18774.0, 18894.0, 18890.0, 18721.0, 19148.0, 18816.0, 19106.0, 
18663.0, 18771.0, 18745.0, 18283.0, 18184.0, 17374.0, 17749.0, 17498.0, 17017.0, 17165.0, 
17141.0, 17216.0, 17240.0, 16993.0, 17037.0, 17304.0, 17295.0, 17172.0, 16884.0, 17105.0, 
17008.0, 16769.0, 16430.0, 16209.0, 15955.0, 16262.0, 15849.0, 15650.0, 15735.0, 15926.0, 
15765.0, 16002.0, 15802.0, 15655.0, 16015.0, 15470.0, 15389.0, 15118.0, 14968.0, 14149.0, 
14508.0, 14446.0, 14129.0, 14312.0, 13939.0, 14260.0, 14206.0, 14144.0, 14160.0, 14274.0, 
14319.0, 14254.0, 14242.0, 13922.0, 13882.0, 13434.0, 13560.0, 13306.0, 13426.0, 13192.0, 
12959.0, 13118.0, 12690.0, 12762.0, 13267.0, 12845.0, 12948.0, 13066.0, 13330.0, 12936.0, 
12780.0, 12705.0, 12899.0, 12737.0, 12363.0, 11828.0, 11873.0, 12123.0, 11868.0, 11705.0, 
11817.0, 11787.0, 11645.0, 11809.0, 11758.0, 11556.0, 11257.0, 11315.0, 11244.0, 10906.0, 
10764.0, 10456.0, 10494.0, 10387.0, 10146.0, 10316.0, 10338.0, 10381.0, 10246.0, 10181.0, 
10098.0, 9732.0, 10120.0, 10096.0, 9610.0, 9654.0, 9428.0, 9320.0, 9355.0, 9163.0, 
9111.0, 8878.0, 8891.0, 8668.0, 8547.0, 8252.0, 8659.0, 8448.0, 8621.0, 8429.0, 
8480.0, 8415.0, 8059.0, 8235.0, 7970.0, 7999.0, 7663.0, 7547.0, 7456.0, 7567.0, 
7640.0, 7617.0, 7633.0, 7425.0, 7339.0, 7552.0, 7323.0, 7383.0, 7274.0, 7292.0, 
7053.0, 6785.0, 6739.0, 6473.0, 6394.0, 6167.0, 6036.0, 6091.0, 5964.0, 5850.0, 
5783.0, 5549.0, 5315.0, 5564.0, 5113.0, 4995.0, 4912.0, 4753.0, 4722.0, 4195.0, 
4246.0, 4077.0, 4066.0, 3821.0, 3819.0, 3521.0, 3498.0, 3264.0, 3296.0, 3194.0, 
3102.0, 2985.0, 3094.0, 2878.0, 2912.0, 2957.0, 2985.0, 2798.0, 2916.0, 2946.0, 
2902.0, 2718.0, 2673.0, 2672.0, 2759.0, 2638.0, 2638.0, 2593.0, 2769.0, 2794.0, 
2883.0, 2965.0, 3124.0, 3063.0, 3218.0, 3049.0, 3191.0, 3189.0, 3176.0, 3206.0, 
3366.0, 3275.0, 3322.0, 3187.0, 3309.0, 3128.0, 3278.0, 3114.0, 3261.0, 3054.0, 
3195.0, 3269.0, 3116.0, 3090.0, 2960.0, 3121.0, 2933.0, 2952.0, 2837.0, 2829.0, 
2885.0, 2704.0, 2510.0, 2369.0, 2414.0, 2464.0, 2369.0, 2297.0, 2341.0, 2248.0, 
2113.0, 2293.0, 2272.0, 2189.0, 2178.0, 2244.0, 2209.0, 2256.0, 2202.0, 2188.0, 
2014.0, 1984.0, 2190.0, 2167.0, 2195.0, 2260.0, 2240.0, 2262.0, 2355.0, 2370.0, 
2263.0, 2418.0, 2378.0, 2413.0, 2335.0, 2483.0, 2270.0, 2295.0, 2151.0, 2262.0, 
2228.0, 2073.0, 2102.0, 2132.0, 2154.0, 1879.0, 1810.0, 1833.0, 1843.0, 1711.0, 
1819.0, 1610.0, 1633.0, 1659.0, 1677.0, 1515.0, 1617.0, 1348.0, 1598.0, 1478.0, 
1428.0, 1416.0, 1502.0, 1397.0, 1353.0, 1453.0, 1424.0, 1423.0, 1334.0, 1347.0, 
1419.0, 1394.0, 1387.0, 1521.0, 1508.0, 1453.0, 1441.0, 1556.0, 1655.0, 1708.0, 
1673.0, 1772.0, 1665.0, 1683.0, 1833.0, 1672.0, 1925.0, 1834.0, 1788.0, 1729.0, 
1625.0, 1728.0, 1653.0, 1623.0, 1608.0, 1442.0, 1343.0, 1454.0, 1425.0, 1254.0, 
1190.0, 1320.0, 1381.0, 1255.0, 1156.0, 1219.0, 1042.0, 943.0, 986.0, 1120.0, 
932.0, 1043.0, 1032.0, 1147.0, 859.0, 923.0, 865.0, 921.0, 801.0, 969.0, 
904.0, 994.0, 958.0, 884.0, 942.0, 958.0, 915.0, 957.0, 867.0, 819.0, 
944.0, 859.0, 830.0, 940.0, 936.0, 871.0, 1039.0, 1022.0, 789.0, 790.0, 
732.0, 861.0, 764.0, 857.0, 761.0, 720.0, 719.0, 805.0, 678.0, 712.0, 
595.0, 489.0, 577.0, 532.0, 467.0, 414.0, 381.0, 315.0, 270.0, 297.0, 
]
end
# A default NVCenter Spectrum Interpolated Curve
const nv_spectrum_spline = SmoothingSplines.fit(SmoothingSplines.SmoothingSpline, NVspectrum.wavelength, NVspectrum.spectrum, 250.)
