'''
exec(open('singleSlopeRunup.py').read())
'''
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as ticker
from scipy.interpolate import LinearNDInterpolator

plt.rc("font", family = "Times New Roman", size = 12.0)
# plt.rcParams['text.usetex'] = True
# plt.rcParams.update({
#     "text.usetex": True,
#     "font.family": "Times New Roman"
# })


class myClassSolitaryWaveInfo:
    def __init__(self, waveHeight=0.1, waterDepth=1.0):
        self.waterDepth = waterDepth
        self.waveHeight = waveHeight

class myInterfaceHeight(myClassSolitaryWaveInfo):
    def __init__(self, inDir, inFile, waveHeight=0.1, waterDepth=1.0):
        myClassSolitaryWaveInfo.__init__(self,waveHeight=waveHeight, waterDepth=waterDepth)

        thisFile = os.path.join(inDir, inFile)
        print(thisFile)
        tmpData = np.loadtxt(
                        open(thisFile, "r"),
                        # delimiter=",",
                        # skiprows=1,
                        # usecols=[0, 2]
                    )
        index = -1
        index += 1; self.time = tmpData[:, index]
        index += 1; self.eta  = tmpData[:, index+1::2]
        self.nSta = np.shape(self.eta)[1]
        self.effectiveH = np.zeros((self.nSta, 2), dtype=float)

        if '00' in inDir:
            pass
        elif '01' in inDir:
            tmpStaGroupIndex = 0
        elif '12' in inDir:
            tmpStaGroupIndex = 1
        elif '23' in inDir:
            tmpStaGroupIndex = 2
        elif '34' in inDir:
            tmpStaGroupIndex = 3
        elif '45' in inDir:
            tmpStaGroupIndex = 4
        
        for i in range(self.nSta):
            self.effectiveH[i, 0] = tmpStaGroupIndex * 10 + i + 0.01
            self.effectiveH[i, 1] = np.max(self.eta[:, i])
        
        self.meanH = np.mean(self.effectiveH[1:-1,1])
        
        
        
waterDepth = 0.5
waveHeight = 0.3 * waterDepth

inFile = 'height.dat'
dirIndex = -1

dirIndex += 1
inDir = f'../postProcessing/interfaceHeight{dirIndex}{dirIndex+1}/0/'
tmpData01 = myInterfaceHeight(inDir, inFile, waterDepth=waterDepth, waveHeight=waveHeight)
np.savetxt(f'effectiveH{dirIndex}{dirIndex+1}.txt', tmpData01.effectiveH, fmt='%1.4e')

dirIndex += 1
inDir = f'../postProcessing/interfaceHeight{dirIndex}{dirIndex+1}/0/'
tmpData12 = myInterfaceHeight(inDir, inFile, waterDepth=waterDepth, waveHeight=waveHeight)
np.savetxt(f'effectiveH{dirIndex}{dirIndex+1}.txt', tmpData12.effectiveH, fmt='%1.4e')

dirIndex += 1
inDir = f'../postProcessing/interfaceHeight{dirIndex}{dirIndex+1}/0/'
tmpData23 = myInterfaceHeight(inDir, inFile, waterDepth=waterDepth, waveHeight=waveHeight)
np.savetxt(f'effectiveH{dirIndex}{dirIndex+1}.txt', tmpData23.effectiveH, fmt='%1.4e')

dirIndex += 1
inDir = f'../postProcessing/interfaceHeight{dirIndex}{dirIndex+1}/0/'
tmpData34 = myInterfaceHeight(inDir, inFile, waterDepth=waterDepth, waveHeight=waveHeight)
np.savetxt(f'effectiveH{dirIndex}{dirIndex+1}.txt', tmpData34.effectiveH, fmt='%1.4e')

dirIndex += 1
inDir = f'../postProcessing/interfaceHeight{dirIndex}{dirIndex+1}/0/'
tmpData45 = myInterfaceHeight(inDir, inFile, waterDepth=waterDepth, waveHeight=waveHeight)
np.savetxt(f'effectiveH{dirIndex}{dirIndex+1}.txt', tmpData45.effectiveH, fmt='%1.4e')

eH = np.mean(np.array([tmpData12.meanH, tmpData23.meanH, tmpData34.meanH])) / waveHeight


dirOut = './'
if True:

    FigureName = os.path.join(dirOut, 'waveProfile')


    FigureNx, FigureNy = 1, 1
    FigureSize = (9, 4)

    # FigureNx, FigureNy = 3, 1
    # FigureSize = (8, 8)

    # FigureNx, FigureNy = 1, 1
    # FigureSize = (8, 4)

    # FigureNx, FigureNy = 2, 2
    # FigureSize = (9, 6)

    # FigureNx, FigureNy = 3, 2
    # FigureSize = (9, 8)

    fig, axs = plt.subplots(ncols=FigureNy, nrows=FigureNx, 
                            figsize=FigureSize,
                            layout="constrained",
                            sharex=False, sharey=False,
                            )
    
    # 
    for FigureIx in range(FigureNx):
        for FigureIy in range(FigureNy):

            FigureIndex = (FigureIx -1 + 1 ) * FigureNy + FigureIy+1
            print(f'ix={FigureIx}, iy={FigureIy}, index={FigureIndex}')
            if FigureNx == FigureNy == 1:
                myAx = axs
            elif FigureNx == 1:
                myAx = axs[FigureIy]
            elif FigureNy == 1:
                myAx = axs[FigureIx]
            else:
                myAx = axs[FigureIx, FigureIy]
            
            myAx.grid()

            textLabelY = "$\eta$/H (-)"
            textLabelX = "Time (s)"  
            # data involking
            if FigureIndex == 0:
                pass

            elif FigureIndex == 1:
                # myAx.set_aspect(1.25)
                textTitle = f'$H$={waveHeight} m, $h$={waterDepth} m' + ", mean.$\\varepsilon=$%.3f"%(eH)
            #  r'\frac{-e^{i\pi}}{2^n}$!'
                if FigureNx * FigureNy > 1:
                    textTitle = f'({FigureIndex}): ' + textTitle
                    
                tmpData = tmpData01
                for iSta in range(0+0, tmpData.nSta, 1):
                    tmp_x = tmpData.time
                    tmp_y = tmpData.eta[:, iSta] / tmpData.waveHeight
                    if iSta == tmpData.nSta-1:
                        myAx.plot(tmp_x, tmp_y, 'k-',
                                label = '$\epsilon$=%.3f'%(tmpData.meanH/tmpData.waveHeight),
                                linewidth=1.5, markersize=12)
                    else:
                        myAx.plot(tmp_x, tmp_y, 'k-',
                                linewidth=1.5, markersize=12)
                
                tmpData = tmpData12
                for iSta in range(0, tmpData.nSta, 1):
                    tmp_x = tmpData.time
                    tmp_y = tmpData.eta[:, iSta] / tmpData.waveHeight
                    if iSta == tmpData.nSta-1:
                        myAx.plot(tmp_x, tmp_y, 'b-',
                                label = '$\epsilon$=%.3f'%(tmpData.meanH/tmpData.waveHeight),
                                linewidth=1.5, markersize=12)
                    else:
                        myAx.plot(tmp_x, tmp_y, 'b-',
                                linewidth=1.5, markersize=12)

                tmpData = tmpData23
                for iSta in range(0, tmpData.nSta, 1):
                    tmp_x = tmpData.time
                    tmp_y = tmpData.eta[:, iSta] / tmpData.waveHeight
                    if iSta == tmpData.nSta-1:
                        myAx.plot(tmp_x, tmp_y, 'm-',
                                label = '$\epsilon$=%.3f'%(tmpData.meanH/tmpData.waveHeight),
                                linewidth=1.5, markersize=12)
                    else:
                        myAx.plot(tmp_x, tmp_y, 'm-',
                                linewidth=1.5, markersize=12)
                
                tmpData = tmpData34
                for iSta in range(0, tmpData.nSta, 1):
                    tmp_x = tmpData.time
                    tmp_y = tmpData.eta[:, iSta] / tmpData.waveHeight
                    if iSta == tmpData.nSta-1:
                        myAx.plot(tmp_x, tmp_y, 'r-',
                                label = '$\epsilon$=%.3f'%(tmpData.meanH/tmpData.waveHeight),
                                linewidth=1.5, markersize=12)
                    else:
                        myAx.plot(tmp_x, tmp_y, 'r-',
                                linewidth=1.5, markersize=12)
                
                tmpData = tmpData45
                for iSta in range(0, tmpData.nSta, 1):
                    tmp_x = tmpData.time
                    tmp_y = tmpData.eta[:, iSta] / tmpData.waveHeight
                    if iSta == tmpData.nSta-2:
                        myAx.plot(tmp_x, tmp_y, 'c-',
                                label = '$\epsilon$=%.3f'%(tmpData.meanH/tmpData.waveHeight),
                                linewidth=1.5, markersize=12)
                    else:
                        myAx.plot(tmp_x, tmp_y, 'c-',
                                linewidth=1.5, markersize=12)
                
             


                if False:
                    FigureXmin, FigureXmax = np.min(tmp_x), np.max(tmp_x)
                    FigureYmin, FigureYmax = np.min(tmp_y), np.max(tmp_y)
                    FigureTickXN = 10
                    FigureTickYN = 10
                    FigureDx = (FigureXmax - FigureXmin) / FigureTickXN
                    FigureDy = (FigureYmax - FigureYmin) / FigureTickYN
                else:
                    FigureXmin, FigureXmax = 0, 20
                    FigureYmin, FigureYmax = -0.2, 0.8
                    FigureYmin, FigureYmax = -0.2, 1.1


                    FigureTickXN = 20
                    FigureTickYN = 13
                    FigureDx = (FigureXmax - FigureXmin) / FigureTickXN
                    FigureDy = (FigureYmax - FigureYmin) / FigureTickYN

                # end
                
            

            # * * * * + * * * * + * * * * + * * * * + * * * * + * * * * + * * * * +
            myAx.legend(
                markerscale = 1.0,
                loc = "lower left",
                fontsize = "small",
                edgecolor = "black",
                facecolor = "white",
                ncol = 5,
                # title = "Velocity CDF",
            )

            myAx.set_title(textTitle, loc = "left")
            '''
                WARN: following lines can not be modified, if YOU know how to do them.
            '''
            # - - - - + - - - - + - - - - + - - - - + - - - - + - - - - + - - - - +
            myAx.set_xlim(FigureXmin, FigureXmax)
            FigureTickX = np.linspace(FigureXmin, FigureXmax, num=FigureTickXN+1)
            # print(f'FigureTickX = {FigureTickX}')
            myAx.set_xticks(FigureTickX)
            FigureTickXLabel = ['%04d'%i for i in FigureTickX] if np.mod(FigureDx,1) == 0 else ['%.1f'%i for i in FigureTickX]
            # myAx.set_xticklabels(FigureTickXLabel)
            if FigureIy == 1 - 1:
                myAx.set_ylabel(textLabelY)
            # end

            myAx.set_ylim(FigureYmin, FigureYmax)
            FigureTickY = np.linspace(FigureYmin, FigureYmax, num=FigureTickYN+1)
            myAx.set_yticks(FigureTickY)
            FigureTickYLabel = ['%04d'%i for i in FigureTickY] if np.mod(FigureDy,1) == 0 else ['%.1f'%i for i in FigureTickY]
            # myAx.set_yticklabels(FigureTickYLabel)
            if FigureIx == FigureNx-1:
                myAx.set_xlabel(textLabelX)
            # end
            # - - - - + - - - - + - - - - + - - - - + - - - - + - - - - + - - - - +
            
            
    # plt.tight_layout()
    # plt.savefig(f'{FigureName}200.png', bbox_inches='tight', dpi=200)
    # plt.savefig(f'{FigureName}400.png', bbox_inches='tight', dpi=400)
    # plt.savefig(f'{FigureName}600.png', bbox_inches='tight', dpi=600)
    plt.savefig(f'{FigureName}.png', bbox_inches='tight', dpi=500)

    # plt.savefig(f'{FigureName}.pdf', bbox_inches='tight')
    plt.show()
# end


# 
