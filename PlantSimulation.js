	// simulateStorageHeatExchange(this.cycleData, false, 0.2, 3600) => execution time < 2s
//firstChargeDischargeSimulation()
//new PlantSimulation(cycleData);
class PlantSimulation {
//function createSimulation(noUI, stepLen, rockHeatingTime) {
	constructor(cycleData, stepLen, rockHeatingTime)
	{
		if (!metrics) return false;
		const swaptHighLpAndtLowHp = metrics.topFillUsedAsColdStorage || metrics.topFillWithIceWithStorage;
	
		this.cycleData = cycleData;
		//this.ignoreSideWalls = getBooleanParam('ignoreSideWalls', true);
		//const hLowLpDischarge = Module.PropsSI('H', 'P', this.cycleData.coldCycle[1].pOut, 'Q', 1, name); 
		this.pControlArea = 0;
		this.pControlVolume = 0;
		this.rPControl = 0;
		this.tHighLp = swaptHighLpAndtLowHp ? this.cycleData.coldCycle[0].tOut : this.cycleData.hotCycle[1].tIn;
		this.tLowHp = swaptHighLpAndtLowHp ?  this.cycleData.hotCycle[1].tIn : this.cycleData.coldCycle[0].tOut;
		this.tHighHp = this.cycleData.hotCycle[1].tOut;
		this.pLowCompressor = (metrics.pLowInCharge ? metrics.pLowInCharge : metrics.pLow);
		this.pHighCompressor = metrics.pHighInCharge ? metrics.pHighInCharge : metrics.pHigh;
		this.tLowLp = Module.PropsSI('T', 'P', this.pLowCompressor, 'Q', 1, name)+0.1;	// 0.01 => crash
		this.pLow = metrics.pLowInDischarge; // metrics.pLow;
		this.pHigh = metrics.pHighInDischarge; // metrics.pLow;
		this.hLowLp = Module.PropsSI('H', 'P', this.pLowCompressor, 'T', this.tLowLp, name);
		this.hHighLp = Module.PropsSI('H', 'P', this.pLowCompressor, 'T', this.tHighLp, name);
		this.hLowHp = Module.PropsSI('H', 'P', this.pHighCompressor, 'T', this.tLowHp, name);
		this.hHighHp = Module.PropsSI('H', 'P', this.pHighCompressor, 'T', this.tHighHp, name);
		this.crushedIceOrRockPorosity = metrics.externalIceStorage ? metrics.crushedIcePorosity : metrics.lpGravelPorosity;
		var idleTime = getFloatParam('idleTime', 50)/100;
		if (idleTime > 0.99  || idleTime < 0) {
			idleTime = idleTime < 0 ? 0 : 0.99;
			$("#idleTime").val(idleTime*100);
		}
		metrics.idleTimeAfter = this.idleTimeAfter = idleTime/(1-idleTime);
		this.hIce = metrics.externalCryoStorage ? this.cycleData.iceStorage.height : 0;
		this.hLpGravel = metrics.hBelowActual;
		this.hFrom = 0; 
		this.gasLpWeightHalf = 0, this.gasHpWeightHalf = 0;
		this.lpTemps = [];
		this.hpTemps = [];
		this.lpHeCCache = {};
		this.lpHeCache = {};
		this.hpHeCCache = {};
		this.hpHeCache = {};
		this.gasMaxWeight = 0;
		this.gasWeight2 = 0;
		this.gasWeightMin = 0;
		this.gasWeightMin2 = 0;
		this.gasLpWeight = 0;
		this.gasWeightMaxDiff = 0;
		this.stepLen = stepLen ? stepLen : getFloatParam("heSimulationthis.stepLen", 0.2);	// the height (m) of each calclulated storace slice
		this.rockHeatingTime = rockHeatingTime ? rockHeatingTime : getIntParam("heSimulationTimeResolution", 900);	// the resolution time in seconds
		this.lpMaxLowDiff = 	getIntParam("lpMaxLowDiff", 50);
		this.lpMaxHighDiff = 	getIntParam("lpMaxHighDiff", getIntParam("maxHighDiff",100));
		this.hpMaxLowDiff = 	getIntParam("hpMaxLowDiff", swaptHighLpAndtLowHp ? 200 : 50);
		this.hpMaxHighDiff = 	getIntParam("hpMaxHighDiff", getIntParam("maxHighDiff",100));		// die fast??
		this.calcHeatExchangeCoeffcient = getBooleanParam("calcHeatExchangeCoeffcient", false );
		this.chargeSnapshots = [];
		this.hpSideIns = {};
		this.lpSideIns = {};
	
		{
			if (Number.isNaN(this.tHighLp) || Number.isNaN(this.tHighHp) || Number.isNaN(this.tLowHp) || Number.isNaN(this.tLowLp)) throw new Error("Debug BP");

			// heat exchange with massive side wall insulations
			const sideRatioMain = 1/Math.sin(Math.PI*metrics.dome.maxAngle/180);	// 90 degrees == 1, 0 => infinity
			// constant 10 simulation heat exchange steps
			this.lpSideIns.width = metrics.insLpSides;
			this.lpSideIns.mainSideWeight = sideRatioMain*metrics.insLpSides*this.stepLen/10*metrics.insGravelDensity*1000;
			this.lpSideIns.uValOfStep = sideRatioMain*this.stepLen*metrics.uOfGravelrockwoolMix/(metrics.insLpSides/10);
			this.lpSideIns.tAmbient = metrics.tAmbient;
			
			this.hpSideIns.width = metrics.ins;
			this.hpSideIns.mainSideWeight = sideRatioMain*metrics.ins*this.stepLen/10*metrics.insGravelDensity*1000;
			this.hpSideIns.uValOfStep = sideRatioMain*this.stepLen*metrics.uOfGravelrockwoolMix/(metrics.ins/10);
			this.hpSideIns.tAmbient = metrics.tAmbient;
//if (Number.isNaN(lpSideIns.mainSideWeight)) throw new Error();
			const sideWallHeatTransfer = function (ctx, tRock, rInner, duration, temps) {
				if (temps.length == 0) {
					const tDiff = (tRock - metrics.tAmbient)/10;
					var tCur = tRock - tDiff/2;
					for (i = 0; i < 10; i++) {
						temps.push(tCur);
						tCur -= tDiff;
//if (!tCur || tCur <= 0 || Number.isNaN(tCur)) throw new Error();
					}
					return 0; // no heat flow in the initialization
				}
				var i = temps.length-1;
				var heatTransferIn, len = 2*Math.PI*(ctx.width + rInner);
				var heatTransferOut = (temps[i] - metrics.tAmbient)*2*ctx.uValOfStep*len*duration;
				//const heatTransferOutWm2 = heatTransferOut /(len*this.stepLen*duration);
				while (i >= 0) {
					const cp = getRockHeatCapasityJinKgK(temps[i]);
					const uVal = i == 0 ? ctx.uValOfStep*2 : ctx.uValOfStep;
					const r = ctx.width*(i + 0.5)/10 + rInner;
					len = 2*Math.PI*r;
					heatTransferIn = ((i == 0 ? tRock : temps[i-1]) - temps[i])*uVal*len*duration;
//if (Number.isNaN(heatTransferIn)) throw new Error("i="+i+", cp="+cp);
					const t = (heatTransferIn-heatTransferOut)/(cp*len*ctx.mainSideWeight);
if (!t || Number.isNaN(t) || typeof heatTransferIn === 'undefined') throw new Error("i="+i);
					temps[i] += t;
					heatTransferOut = heatTransferIn;
					i--;
				}
/*	
				This logging shows now the expected results: in stable state heat flow in and out of insulation is 20 W/m2,
				and max out heat flow is the config param 20 W/m2
				if (heatTransferIn > maxHeatTransferIn || heatTransferIn < minHeatTransferIn) {
					if (heatTransferIn > maxHeatTransferIn) maxHeatTransferIn = heatTransferIn;
					else minHeatTransferIn = heatTransferIn;
					console.log("tRock="+tRock+", in W/m2="+(heatTransferIn/(len*this.stepLen*duration))+", out W/m2="+heatTransferOutWm2+", temps="+JSON.stringify(temps));
				}
*/
if (typeof heatTransferIn === 'undefined' || Number.isNaN(heatTransferIn)) throw new Error("i="+i);
				return heatTransferIn; 	// heat transferred in Joules
			}
console.log("lpSideIns="+JSON.stringify(this.lpSideIns));
console.log("this.hpSideIns="+JSON.stringify(this.hpSideIns));
			this.hpSideIns.sideWallHeatTransfer = sideWallHeatTransfer;
			this.lpSideIns.sideWallHeatTransfer = sideWallHeatTransfer;
			
			// external cryo storage exists => create heat exchange for its massive circular side wall 
			if (this.hIce) {
				// ice side wall:  inner gravel+stone wool insulation + land fill + external insulation (stone wool + XPS)
				const uForzenLandFill = getFloatParam("uForzenLandFill", 3.0);
				const frozenLandFillDensity = getFloatParam("frozenLandFillDensity", 2.1);
				const cInnerIns = 5;
				const cLandFills = 5;	// ~12 massive land fills
				const sideRatioIce = 1/Math.sin(Math.PI/4);	// constant 45 degrees

				// HeatFlowOut = uOfGravelrockwoolMix*deltaT/len =>
				const innerCryoIns = metrics.uOfGravelrockwoolMix * (metrics.tAmbient - this.tLowLp) / metrics.maxHeatFlowOut; // width of inner insulataion
				const extInnerIceInsWeight = sideRatioIce*innerCryoIns*this.stepLen/cInnerIns*metrics.insGravelDensity*1000;
				const uValOfInnerIceStep = sideRatioIce*this.stepLen*metrics.uOfGravelrockwoolMix/(innerCryoIns/cInnerIns);
				const landFillWidth = (this.cycleData.iceStorage.sideWidth/sideRatioIce - innerCryoIns);
				const landFillInsWeight = sideRatioIce*landFillWidth*this.stepLen/cLandFills*frozenLandFillDensity*1000;
				const uValOfIceLandFillStep = sideRatioIce*this.stepLen*uForzenLandFill/(landFillWidth/cLandFills);
				const uValOfIceLastStep = sideRatioIce*this.stepLen*(metrics.tAmbient - this.tLowLp) / metrics.maxHeatFlowOut + uValOfIceLandFillStep*2;
				const uTotal = uValOfInnerIceStep*cInnerIns + uValOfIceLandFillStep*cLandFills + uValOfIceLastStep;
				
				const extSideWallHeatTransfer = function (ctx, tRock, rInner, duration, temps) {
					if (temps.length == 0) {
						const tDiff = (metrics.tAmbient - tRock);
						var tCur = tDiff*uValOfIceLastStep/uTotal;
						const tDiffLandFill = tDiff*uValOfIceLandFillStep/uTotal;
						const tDiffInner = tDiff*uValOfInnerIceStep/uTotal;
						temps.push(tCur);
						temps.push(tCur += tDiffLandFill);
						temps.push(tCur += tDiffLandFill);
						temps.push(tCur += tDiffLandFill);
						temps.push(tCur += tDiffLandFill);
						temps.push(tCur += (tDiffLandFill + tDiffInner)/2);
						temps.push(tCur += tDiffInner);
						temps.push(tCur += tDiffInner);
						temps.push(tCur += tDiffInner);
						temps.push(tCur += tDiffInner);
						temp = temp.reverse();	
						return 0; // no heat flow in the initialization
					}
					var i = temps.length-1;
					var heatTransferIn, len = 2*Math.PI*(landFillWidth + innerCryoIns + rInner);
					var heatTransferOut = (temps[i] - metrics.tAmbient)*uValOfIceLastStep*len*duration;
					const heatTransferOutWm2 = heatTransferOut /(len*this.stepLen*duration);
					while (i >= 0) {
						const cp = getRockHeatCapasityJinKgK(temps[i]);		// frozen landfill probably have about the same cp as rock in cryo temps
						const r = i >= 5 ? landFillWidth*(i - 4.5)/5 + innerCryoIns : innerCryoIns * (i+0.5)/5 + rInner;
						len = 2*Math.PI*r;
						const uVal = i > 5 ? uValOfIceLandFillStep : 
									 i == 5 ? 1/(1/(2*uValOfIceLandFillStep)+1/(2*uValOfInnerIceStep)) :	// 1/k = 1/k1 + 1/k2
									 i > 0 ? uValOfInnerIceStep : uValOfInnerIceStep*2;
						const weight = i >= 5 ? landFillInsWeight : extInnerIceInsWeight;
						heatTransferIn = ((i == 0 ? tRock : temps[i-1]) - temps[i])*uVal*len*duration;
						const t = (heatTransferIn-heatTransferOut)/(cp*len*weight);
						temps[i] += t;
if (Number.isNaN(heatTransferIn) || Number.isNaN(t)) throw new Error();
						heatTransferOut = heatTransferIn;
						i--;
					}
if (typeof heatTransferIn === 'undefined' || Number.isNaN(heatTransferIn)) throw new Error("i="+i);
					return heatTransferIn;	// heat transferred in Joules
				}
				this.lpSideIns.extSideWallHeatTransfer = extSideWallHeatTransfer;
			}
		}
		
		if (!this.lpHeCCache.size) {
			this.lpHeCCache.stepLen = this.stepLen;
			//this.lpHeCCache.liquidStorage = lpLiquidStorage;
			this.lpHeCCache.size = (Math.floor(this.tHighLp) - Math.floor(this.tLowLp) + 10);
			this.lpHeCCache.tMin = Math.floor(this.tLowLp);
			this.lpHeCCache.tLow = this.tLowLp;
			this.lpHeCCache.tHigh = this.tHighLp;
			this.lpHeCCache.hLow = this.hLowLp;
			this.lpHeCCache.hHigh = this.hHighLp;
			this.lpHeCCache.p = this.pLowCompressor;
			this.lpHeCCache.pLowDischarge = this.pLow;
			this.lpHeCCache.tSat = Module.PropsSI('T', 'P', this.lpHeCCache.p, 'Q', 1, name);
			this.lpHeCCache.hSat = Module.PropsSI('H', 'P', this.lpHeCCache.p, 'Q', 1, name);
console.log("lpHeCCache: tSat="+this.lpHeCCache.tSat+", hSat="+this.lpHeCCache.hSat+", tLow="+this.lpHeCCache.tLow+", hLow="+this.lpHeCCache.hLow+", p="+this.lpHeCCache.p);
			this.lpHeCCache.porosity = metrics.lpGravelPorosity;
			this.lpHeCCache.gravelMin = metrics.gravelMinLp / 1000;
			this.lpHeCCache.gravelMax = metrics.gravelMaxLp / 1000;
			this.lpHeCCache.medianDP = (this.lpHeCCache.gravelMin + this.lpHeCCache.gravelMax) / 2;
//console.log("lpHeCCache.gravelMin="+this.lpHeCCache.gravelMin+", metrics.gravelMaxLp="+metrics.gravelMaxLp+", this.lpHeCCache.medianDP="+this.lpHeCCache.medianDP);
			this.lpHeCCache.gravelDensityKgInM3 = metrics.lpGravelDensity * 1000;
			this.lpHeCCache.rockHeatingTime = this.rockHeatingTime;	// recalc rock temperature in 1 h steps
			this.lpHeCCache.sideIns = this.lpSideIns;
			this.lpHeCCache.iceLength = this.hIce/this.stepLen;
		}
		if (!this.lpHeCache.size) {
			this.lpHeCache.stepLen = this.stepLen;
			//this.lpHeCCache.liquidStorage = lpLiquidStorage;
			this.lpHeCache.size = (Math.floor(this.tHighLp) - Math.floor(this.tLowLp) + 10);
			this.lpHeCache.tMin = Math.floor(this.tLowLp);
			this.lpHeCache.tLow = this.tLowLp;
			this.lpHeCache.tHigh = this.tHighLp;
			this.lpHeCache.hLow = this.hLowLp;
			this.lpHeCache.hHigh = Module.PropsSI('H', 'P', this.pLow, 'T', this.tHighLp, name);;
			this.lpHeCache.p = this.pLow; 
			this.lpHeCache.pLowCharge = this.pLowCompressor;  
			this.lpHeCache.tSat = Module.PropsSI('T', 'P', this.lpHeCache.p, 'Q', 1, name);
			this.lpHeCache.hSat = Module.PropsSI('H', 'P', this.lpHeCache.p, 'Q', 1, name);
console.log("lpHeCache: tSat="+this.lpHeCache.tSat+", hSat="+this.lpHeCache.hSat+", tLow="+this.lpHeCache.tLow+", hLow="+this.lpHeCache.hLow+", p="+this.lpHeCache.p);
			this.lpHeCache.porosity = metrics.lpGravelPorosity;
			this.lpHeCache.gravelMin = metrics.gravelMinLp / 1000;
			this.lpHeCache.gravelMax = metrics.gravelMaxLp / 1000;
			this.lpHeCache.medianDP = (this.lpHeCache.gravelMin + this.lpHeCache.gravelMax) / 2;
//console.log("lpHeCache.gravelMin="+this.lpHeCache.gravelMin+", metrics.gravelMaxLp="+metrics.gravelMaxLp+", this.lpHeCache.medianDP="+this.lpHeCache.medianDP);
			this.lpHeCache.gravelDensityKgInM3 = metrics.lpGravelDensity * 1000;
			this.lpHeCache.rockHeatingTime = this.rockHeatingTime;	// recalc rock temperature in 1 h steps
			this.lpHeCache.sideIns = this.lpSideIns;
			this.lpHeCCache.iceLength = this.hIce/this.stepLen;
		}
		if (!this.hpHeCache.size) {
			this.hpHeCache.stepLen = this.stepLen;
			//this.lpHeCCache.liquidStorage = hpLiquidStorage;
			this.hpHeCache.size = (Math.floor(this.tHighHp+2) - Math.floor(this.tLowHp) + 10);
			this.hpHeCache.tMin = Math.floor(this.tLowHp) - 1;
			this.hpHeCache.tLow = this.tLowHp;
			this.hpHeCache.tHigh = this.tHighHp;
			this.hpHeCache.hLow = Module.PropsSI('H', 'P', this.pHigh, 'T', this.tLowHp, name);
			this.hpHeCache.hHigh = Module.PropsSI('H', 'P', this.pHigh, 'T', this.tHighHp, name);
			this.hpHeCache.porosity = metrics.hpGravelPorosity;
			this.hpHeCache.p = this.pHigh;
			this.hpHeCache.pHighCharge = this.pHighCompressor;
			this.hpHeCache.gravelMin = metrics.gravelMinHp / 1000;
			this.hpHeCache.gravelMax = metrics.gravelMaxHp / 1000;
			this.hpHeCache.medianDP = (this.hpHeCache.gravelMin + this.hpHeCache.gravelMax) / 2;
			this.hpHeCache.gravelDensityKgInM3 = metrics.gravelDensity * 1000;
			this.hpHeCache.rockHeatingTime = this.rockHeatingTime;	// recalc rock temperature in 1 h steps
			this.hpHeCache.sideIns = this.hpSideIns;
		}
		if (!this.hpHeCCache.size) {
			this.hpHeCCache.stepLen = this.stepLen;
			//this.lpHeCCache.liquidStorage = hpLiquidStorage;
			this.hpHeCCache.size = (Math.floor(this.tHighHp+2) - Math.floor(this.tLowHp) + 10);
			this.hpHeCCache.tMin = Math.floor(this.tLowHp) - 1;
			this.hpHeCCache.tLow = this.tLowHp;
			this.hpHeCCache.tHigh = this.tHighHp;
			this.hpHeCCache.hLow = this.hLowHp;
			this.hpHeCCache.hHigh = this.hHighHp;
			this.hpHeCCache.porosity = metrics.hpGravelPorosity;
			this.hpHeCCache.p = this.pHighCompressor;
			this.hpHeCCache.pHighDischarge = this.pHigh;
			this.hpHeCCache.gravelMin = metrics.gravelMinHp / 1000;
			this.hpHeCCache.gravelMax = metrics.gravelMaxHp / 1000;
			this.hpHeCCache.medianDP = (this.hpHeCCache.gravelMin + this.hpHeCCache.gravelMax) / 2;
			this.hpHeCCache.gravelDensityKgInM3 = metrics.gravelDensity * 1000;
			this.hpHeCCache.rockHeatingTime = this.rockHeatingTime;	// recalc rock temperature in 1 h steps
			this.hpHeCCache.sideIns = this.hpSideIns;
		}

		this.getInitialRockTemp = function(h, hMid, hMax, heLength, tMin, tMax ) {
			const heTop = hMid + heLength/2;
			const heBottom = hMid - heLength/2;
			if (h >= heTop) {
				return tMax; // - 5*(hMax-h)/(hMax-heTop);
			}
			else if (h <= heBottom) {
				return tMin; // + 5*h/heBottom;
			}
			else if (h <= hMid) {
				return tMin + Math.pow((h-heBottom)/(heLength/2),2) * (tMax-tMin)/2;
			}
			else  {
				return tMax - Math.pow((heTop-h)/(heLength/2),2) * (tMax-tMin)/2;
			}
		}
		// metrics.heLengthHp, metrics.heLengthLp,  metrics.hColdHpGravel, metrics.hMidLpFromGround;
		//const iceArea = 1.2*metrics.lowPressureStorageVolume / this.hIce;	// start with something
		
		if (this.lpTemps.length == 0) {
			this.lpInitialRockTemp = [];
			const hMidLp = this.hLpGravel+metrics.hMidLpFromGround;
			var volume = 0;
			//const includeHeatLeaks = getBooleanParam('includeHeatLeaks',true);
			//var isLpBottom = includeHeatLeaks;
			var hLP = metrics.topFillWithIceWithStorage ? this.cycleData.iceStorage.height : this.hIce+this.hLpGravel+(metrics.topFillUsedAsWeight ||metrics.noTopFill ? 0 : metrics.hAbove);
			var iStep = 0;
			for (h = this.stepLen; h < hLP; h += this.stepLen) {
				var item = {};
				// Set lp storage to the min energy level
				//item.tRock = this.getInitialRockTemp(h, hLP/5, hLP, hLP*2/5, lpLiquidStorage.tHigh, this.tHighLp);
				item.tRock =  metrics.topFillWithIceWithStorage ? this.tHighLp : this.tLowLp; // this.getInitialRockTemp(h, hLP*0.1, hLP, metrics.topFillWithIceWithStorage ? 1 : metrics.heLengthLp, this.tLowLp, this.tHighLp);
				//item.tRock = 283; // 10 C
				//item.tRock = h < hMidLp ? this.tLowLp : this.tHighLp;
				//item.hGasOut = h < hMidLp ? this.hLowLp : this.hHighLp;
				
				this.lpInitialRockTemp.push(item.tRock);
				if (metrics.topFillWithIceWithStorage) {
					item.area = this.cycleData.iceStorage.areaM2;
				}
				else if (h < this.hIce) {
					item.area = this.cycleData.iceStorage.getArea(h) - metrics.hasLiquidStorage ? 0 : this.pControlArea;
					item.type = 'ice';
					item.r = Math.sqrt(item.area/Math.PI);
					item.sideWallTemps = [];
					this.lpSideIns.extSideWallHeatTransfer(this.lpSideIns, item.tRock, item.r, this.rockHeatingTime, item.sideWallTemps);
				}
				else {
					var r;
					if (h > this.hLpGravel) {
						const hAboveCur = h - this.hLpGravel;
						r = metrics.rGround - (metrics.rGround-metrics.rTop)*hAboveCur/metrics.hAbove;
						if (!r) {
							console.log(r+" = "+metrics.rGround+" - ("+metrics.rGround+"-"+metrics.rTop+")*"+hAboveCur+"/"+metrics.hAbove);
						}
					}
					else {
						if (false && isLpBottom) {
							logHeatLeak = true;
							// double insulatio in the bottom because of 
							console.log("Heat leak from HP storage: "+Math.round(item.heatLeakW)+" W");
							isLpBottom = false;
						}
						//r = metrics.rOpenBottomActual + (metrics.rGroundBelowActual-metrics.rOpenBottomActual)*(h+metrics.hPressureControl)/metrics.hBelowActual;
						r = metrics.rOpenBottomActual + (metrics.rGroundBelowActual-metrics.rOpenBottomActual)*(h)/metrics.hBelowActual;
					}
					item.area = (r*r*Math.PI);		// correction factor for low heat capasity
					item.r = r;
					item.sideWallTemps = [];
					this.lpSideIns.sideWallHeatTransfer(this.lpSideIns, item.tRock, item.r, this.rockHeatingTime, item.sideWallTemps);
				}
				volume += item.area * this.stepLen;
if (!item.area) {
// Bad area=NaN, r=NaN, h=59.5, metrics.rOpenBottomActual=148.37160254662282, metrics.rGroundBelowActual=166.47297467271602, metrics.hBelowActual=59.2069204167878
	console.log("Bad area="+item.area+", r="+r+", h="+h+", metrics.rOpenBottomActual="+metrics.rOpenBottomActual+", metrics.rGroundBelowActual="+metrics.rGroundBelowActual+", metrics.hBelowActual="+metrics.hBelowActual);
	return false;
}
				this.lpTemps.push(item);
				iStep++;
			}
/*
			const lpCpToGasCp = getAverageRockHeatCapasitykJinKg(this.tLowLp, this.tHighLp)*1000*(this.tHighLp-this.tLowLp)*volume
					/ (Module.PropsSI('H', 'P', pLow, 'T', this.tHighLp, name)-Module.PropsSI('H', 'P', pLow, 'T', this.tLowLp, name));
			const cpHighGravel = getRockHeatCapasityJinKgK(this.tHighLp), cpLowGravel = getRockHeatCapasityJinKgK(this.tLowLp);
			const cpLowIce = getRockHeatCapasityJinKgK(this.tLowLp, 'ice');
			// iceVolume*metrics.lpGravelDensity*cpLowIce+ cpLowGravel*volume*metrics.lpGravelDensity=cpHighGravel*volume*metrics.lpGravelDensity =>
			const iceVolume = (cpHighGravel*volume-cpLowGravel*volume) / cpLowIce;
console.log("LP Storage size:"+myRound(volume*metrics.lpGravelDensity/1000000, 1)+" Mt, lpCpToGasCp="+lpCpToGasCp+", iceVolume="+(iceVolume/1000000));
*/
		}
		else {
			for (i = 0; i < this.lpTemps.length; i++) {
				this.lpInitialRockTemp[i] = this.lpTemps[i].tRock;
			}
		}
//console.log("lpTemps:"+JSON.stringify(this.lpTemps));
		if (this.hpTemps.length == 0) {
			const hMidHp = metrics.hColdHpGravel;
			this.hpInitialRockTemp = [];
			volume = 0;
			//const this.tLowHp= Module.PropsSI('T', 'P', this.hpHeCache.p, 'H', this.hLowHp, name);
			//const this.tHighHp= Module.PropsSI('T', 'P', this.hpHeCache.p, 'H', this.hHighHp, name);
			for (h = this.stepLen; h < metrics.dome.hBottomI+this.stepLen; h += this.stepLen) {
				var item = {};
				//item.tRock = this.getInitialRockTemp(h, metrics.dome.hBottomI-metrics.heLengthHp/2, metrics.dome.hBottomI, metrics.heLengthHp, this.tLowHp, this.tHighHp-this.hpMaxHighDiff);
				item.tRock = metrics.topFillWithIceWithStorage || h < (metrics.dome.hBottomI/2) ? this.tLowHp : this.tHighHp; // this.getInitialRockTemp(h, metrics.dome.hBottomI*0.9, metrics.dome.hBottomI, metrics.heLengthHp, this.tLowHp, this.tHighHp);
				//item.tRock = 283; // 10 C
				//item.tRock = h < hMidHp ? this.tLowHp : this.tHighHp;
				//item.hGasOut = h < hMidHp ? this.hLowHp : this.hHighHp;
				this.hpInitialRockTemp.push(item.tRock);
				const r = metrics.dome.rBottomI + (metrics.dome.rTopI-metrics.dome.rBottomI)*h/metrics.dome.hBottomI;
				item.area = (r*r*Math.PI);
				volume += item.area * this.stepLen;
				item.r = r;
				item.sideWallTemps = [];
				this.hpSideIns.sideWallHeatTransfer(this.hpSideIns, item.tRock, item.r, this.rockHeatingTime, item.sideWallTemps);
				
if (!item.area) {
	console.log("Bad area="+item.area+", r="+r+", h="+h);
	return false;
}
				this.hpTemps.push(item);
			}
			const hpCpToGasCp = getAverageRockHeatCapasitykJinKg(this.tLowHp, this.tHighHp)*1000*(this.tHighHp-this.tLowHp)*volume 
					/ (Module.PropsSI('H', 'P', metrics.pHigh, 'T', this.tHighHp, name)-Module.PropsSI('H', 'P', metrics.pHigh, 'T', this.tLowHp, name));
console.log("HP Storage size:"+myRound(volume*metrics.gravelDensity/1000000, 1)+" Mt, hpCpToGasCp="+hpCpToGasCp);
			//console.log("hpTemps="+JSON.stringify(this.hpTemps));
		}
		else {
			for (i = 0; i < this.hpTemps.length; i++) {
				this.hpInitialRockTemp[i] = this.hpTemps[i].tRock;
			}
		}
		this.tMaxIceStorage = 0;
		this.tMinIceStorage = 0;

		// TBD: randomBalancedDays
//console.log("hpTemps:"+JSON.stringify(this.hpTemps));
		this.simulateCharge = function(calcHC, hours )
		{
			var loop = -1, hour, tGas, height, hGas;
			for (hour = 0; hour < hours; hour++) {
				tGas = this.tLowLp;
				hGas = this.hLowLp;
				//lpLiquidStorage.useStorage(this.cycleData.hotCycle[1].massFlow*1000, hGas, tGas, this.rockHeatingTime, true);
				//tGas = lpLiquidStorage.tOut;
				hGas = 0;
				for (height = this.hFrom; height < this.lpTemps.length; height++) {
					if (!processGasToGravelHeatExchange( this.lpHeCCache, this.lpTemps[height], tGas, hGas, this.cycleData.hotCycle[1].massFlow, calcHC, loop == 0 && height<10 )) {
						console.log(height+":this.lpTemps["+(height-10)+"-"+(height+10)+"]="+JSON.stringify(this.lpTemps.slice(height >= 10 ? height-10 : 0,(height+10) < this.lpTemps.length ? height+10 : this.lpTemps.length)));
						return fromHere(loop,height);
					}
					tGas = this.lpTemps[height].tRock;
					hGas = this.lpTemps[height].hGasOut;
				}
				processGravelHeatTransfers(this.lpHeCCache, this.lpTemps, this.stepLen, (1 + this.idleTimeAfter) * this.lpHeCCache.rockHeatingTime, loop == 0 && (height < 5 || height > this.lpTemps.length-5));
				
				tGas = this.tHighHp;
				hGas = this.hHighHp;
				for (height = this.hpTemps.length-1; height >= 0; height--) {
					if (!processGasToGravelHeatExchange( this.hpHeCCache, this.hpTemps[height], tGas, hGas, this.cycleData.hotCycle[1].massFlow, calcHC, loop == 0 && height > this.hpTemps.length-10 )) {
						console.log(height+":"+JSON.stringify(this.hpTemps));
						return fromHere(loop,height);
					}
					tGas = this.hpTemps[height].tRock;
					hGas = this.hpTemps[height].hGasOut;
					/*if (!tGas || !hGas) {
						console.log(height+": bad item:"+JSON.stringify(this.lpTemps[height]));
						return false;
					}*/
				}
				processGravelHeatTransfers(this.hpHeCCache, this.hpTemps, this.stepLen, (1 + this.idleTimeAfter) * this.hpHeCCache.rockHeatingTime, loop == 0 && (height < 5 || height > this.hpTemps.length-5));

				//
				//	Take the first snap shot of lp and hp gas weights in half charged storage
				// Needed to resolve the floating charge pressure of mid charged storage
				//
				if (this.midHour == hour) {
					var k, d, tPrev = 0, h = this.stepLen/2
					for (k = this.hFrom; k < this.lpTemps.length; k++) {
						if (Math.abs(this.lpTemps[k].tRock - tPrev)>1) {
							d = Module.PropsSI('D', 'P', this.lpHeCCache.p, 'T', this.lpTemps[k].tRock, name);
							tPrev = this.lpTemps[k].tRock;
						}
						if (h < this.hIce) {
							this.gasLpWeightHalf += d*this.lpTemps[k].area*this.stepLen*this.crushedIceOrRockPorosity;
						}
						else {
							this.gasLpWeightHalf += d*this.lpTemps[k].area*this.stepLen*metrics.lpGravelPorosity;
						}
					}
					tPrev = 0;
					for (k = 0; k < this.hpTemps.length; k++) {
						if (Math.abs(this.hpTemps[k].tRock - tPrev)>1)
						{
							d = Module.PropsSI('D', 'P', this.hpHeCCache.p, 'T', this.hpTemps[k].tRock, name);
							tPrev = this.hpTemps[k].tRock;
						}
						this.gasHpWeightHalf += d*this.hpTemps[k].area*this.stepLen*metrics.hpGravelPorosity;
					}
					this.midHour = -1;
				}
				//hpLiquidStorage.useStorage(this.cycleData.hotCycle[1].massFlow*1000, hGas, tGas, this.rockHeatingTime, false);
			}
			this.tMinIceStorage = this.lpTemps[this.hFrom].tRock;
			return true;
		}
		this.simulateDischarge = function(calcHC, hours )
		{
			var loop = -1, hour, tGas, height, hGas;
			for (hour = 0; hour < hours; hour++) {
				tGas = this.tHighLp;
				hGas = this.lpHeCache.hHigh;
				for (height = this.lpTemps.length-1; height >= this.hFrom; height--) {
					if (!processGasToGravelHeatExchange( this.lpHeCache, this.lpTemps[height], tGas, hGas, this.cycleData.hotCycle[0].massFlow, calcHC, loop == 0 && height > this.lpTemps.length-10 )) {
						console.log(height+":this.lpTemps["+(height-10)+"-"+(height+10)+"]="+JSON.stringify(this.lpTemps.slice(height-10,height+10)));
						return fromHere(hour,height*this.stepLen);
					}
					tGas = this.lpTemps[height].tRock;
					hGas = this.lpTemps[height].hGasOut;
				}
				//lpLiquidStorage.useStorage(this.cycleData.hotCycle[0].massFlow*1000, hGas, tGas, this.rockHeatingTime, false, this.lpHeCache.p, this.lpHeCache.tLow);
				processGravelHeatTransfers(this.lpHeCache, this.lpTemps, this.stepLen, this.lpHeCache.rockHeatingTime, loop == 0 && (height < 5 || height > this.lpTemps.length-5));
				tGas = this.tLowHp;
				hGas = this.lpHeCache.hLow;
				//hpLiquidStorage.useStorage(this.cycleData.hotCycle[0].massFlow*1000, hGas, tGas, this.rockHeatingTime, true);
				//tGas = hpLiquidStorage.tOut;
				hGas = 0; // hpLiquidStorage.hOut;
				for (height = 0; height < this.hpTemps.length; height++) {
					if (!processGasToGravelHeatExchange( this.hpHeCache, this.hpTemps[height], tGas, hGas, this.cycleData.hotCycle[0].massFlow, calcHC, loop == 0 && height<10 )) {
						return fromHere(loop,height*this.stepLen);
					}
					tGas = this.hpTemps[height].tRock;
					hGas = this.hpTemps[height].hGasOut;
				}
				processGravelHeatTransfers(this.hpHeCache, this.hpTemps, this.stepLen, this.lpHeCCache.rockHeatingTime, loop == 0 && (height < 5 || height > this.hpTemps.length-5));
			}
			if (metrics.externalCryoStorage) {
				const h = Math.round(this.cycleData.iceStorage.height/this.stepLen);
				if (this.lpTemps[h].tRock > this.tMaxIceStorage) this.tMaxIceStorage = this.lpTemps[h].tRock;
				this.tMax = this.lpTemps[this.lpTemps.length-1].tRock;
			}
			return true;
		}
//console.log("lpTemps:"+JSON.stringify(this.lpTemps, null, 2));
//console.log("hpTemps:"+JSON.stringify(this.hpTemps, null, 2));
		this.drawHeatExchangeDiagram2 = function( id, items, title, bottomLines, arraysOfVals, minMaxDefault, fpCustomizer ) 
		{
			//console.log(id+":"+title+", bottomLines="+JSON.stringify(bottomLines)); //+", arraysOfVals="+JSON.stringify(arraysOfVals));
			var xVals = [];
			var yTitles =  [];
			yTitles.push("K");
			for (i = 1; i < arraysOfVals.length; i++) yTitles.push("");
			for (h = 0; h < items.length; h++) {
				xVals.push((h+0.5)*this.stepLen);
			}
			drawDiagram4Y( id, xVals, 'm', arraysOfVals, yTitles, title, bottomLines, minMaxDefault, fpCustomizer );
		}
		this.drawStoneAndIceEnthalpies = function( id, tLowLp, tHighHp, fpCustomizer) {
			tLowLp = Math.round(tLowLp), tHighHp = Math.round(tHighHp);
			var xVals = [], arraysOfVals = [], yTitles = [], bottomLines = [];
			minMaxDefault = {min : Math.round(tLowLp/10)*10, max : (1+Math.round(tHighHp/10))*10};
			var iceCp = [], rockCp = [], rockH = [], iceH = [], h = 0;
			bottomLines.push('Gravel heat capasity kJ/m3/K');
			bottomLines.push('Crushed ice heat capasity kJ/m3/K');
			bottomLines.push('Gravel enthalpy MJ/m3 ('+Math.round(1000*metrics.lpGravelDensity)+' kg/m3)');
			bottomLines.push('Crushed ice enthalpy MJ/m3 ('+Math.round(1000*metrics.crushedIceDensity)+' kg/m3)');
			const dLiq = Module.PropsSI('D', 'P', this.lpHeCache.p, 'Q', 0, name);
			metrics.liqVolM3 = Math.round((this.gasMaxWeight-this.gasWeightMin)/dLiq);
			if (this.chargeSnapshots) {
				const first = this.chargeSnapshots[0];
				const last = this.chargeSnapshots[this.chargeSnapshots.length-1];
				this.gasWeightMaxDiff = first.curKapasityKg + first.curChangeKg;
				this.gasWeightMinDiff = last.curKapasityKg + last.curChangeKg;
			}
			const extraTitle = "Gas weight in tons min/max: "+Math.round(this.gasWeightMin/1000)+' / '+Math.round(this.gasMaxWeight/1000)+" => Max liquid volume: "+metrics.liqVolM3+" m3";
			const extraTitle2 = "Charge-to-discharge switch time from "+
					Math.round(this.gasWeightMaxDiff/1000/coldCycleG[0].massFlow/60)+" to "+Math.round(this.gasWeightMinDiff/1000/coldCycleG[0].massFlow/60) + " minutes";
			const extraTitle3 = ''; /*"Discharge high/low pressures (kPa): "+Math.round(pHighDischargeActual/1000)+'/'+myRound(pLowDischargeActual/1000)+" - "+
									Math.round(pHighDischargeActual2/1000)+'/'+myRound(pLowDischargeActual2/1000);*/

			yTitles.push('kJ/m3/K');
			yTitles.push('kJ/m3/K');
			yTitles.push('MJ/m3');
			yTitles.push('MJ/m3');
			for (t = this.tLowLp; t <= 273; t++) {
				const cp = getIceHeatCapasityKJinKgK(t)*1000*metrics.crushedIceDensity;	// ice density
				iceCp.push(cp);
				h = h + (cp/1000);
				iceH.push(h);
			}
//console.log(JSON.stringify(iceH));
			h = 0;
			for (t = this.tLowLp; t <= this.tHighHp; t++) {
				xVals.push(t);
				const cpRock = getRockHeatCapasityJinKgK(t) * metrics.lpGravelDensity;
				rockCp.push(cpRock);
				h = h + (cpRock/1000);
				rockH.push(h);
			}
//console.log(JSON.stringify(rockH));
//console.log(JSON.stringify(iceCp));
//console.log(JSON.stringify(rockCp));
//console.log("Ice enthalpy(t): "+Math.round(this.tLowLp+50)+"="+iceH[50]+", "+Math.round(this.tLowLp+100)+"="+iceH[100]+", "+Math.round(this.tLowLp+150)+"="+iceH[150]+", 273="+iceH[Math.round(273-this.tLowLp)]);
//console.log("RockEnthalpy(t): "+Math.round(this.tLowLp+50)+"="+rockH[50]+", "+Math.round(this.tLowLp+100)+"="+rockH[100]+", "+Math.round(this.tLowLp+150)+"="+rockH[150]+", 273="+rockH[Math.round(273-this.tLowLp)]);
			arraysOfVals.push(rockCp);
			arraysOfVals.push(iceCp);
			arraysOfVals.push(rockH);
			arraysOfVals.push(iceH);

			showIceEndEnthalpy = function (ctx, dg) {
				const i = iceH.length-1;
				ctx.font = "12px Arial";
				ctx.fillStyle = 'blue';
				ctx.fillText( Math.round(iceH[i]/metrics.crushedIceDensity)+" kJ/kg", toX(dg, xVals[i]), toY(dg, iceH[i]));
			};
			drawDiagram4Y( id, xVals, 'K', arraysOfVals, yTitles, "Crushed ice and gravel heat capasities and enthalpies in this simulation,"+extraTitle+","+extraTitle2+","+extraTitle3, bottomLines, null, showIceEndEnthalpy );
		}
		this.drawHeatExchangeDiagram = function( id, items, title, initialTemps, minMaxDefault ) 
		{
			var xVals = [], bottomLines = [], arrayOfyVals = [], yTitles =  [], i;

			arrayOfyVals.push([]);
			arrayOfyVals.push([]);
			arrayOfyVals.push([]);
			arrayOfyVals.push([]);
			xVals = [];

			bottomLines.push("Gravel temp after "+chargeDays+" charge days");
			yTitles.push("K");
			bottomLines.push(+chargeDischargeDays+" days 16h charge+8h discharge");
			yTitles.push("");
			bottomLines.push(dischargeDays+" days 16h discharge+8 charge");
			yTitles.push("");
			bottomLines.push("Gravel temp after "+idleDays+" idle days");
			yTitles.push("");

			for (h = 0; h < items.length; h++) {
				xVals.push((h+0.5)*this.stepLen);
				try {
					//arrayOfyVals[0].push(items[h].heCofficientInM3);
					//arrayOfyVals[1].push(items[h].heatTransferCoefficient);
					arrayOfyVals[0].push(items[h].tC);
					arrayOfyVals[1].push(items[h].tCDC);
					arrayOfyVals[2].push(items[h].tDC);
					arrayOfyVals[3].push(items[h].tRock);
				}
				catch (e) {
					console.log("Ignored h="+h+":"+JSON.stringify(items[h]));
				}
			}
	//console.log(JSON.stringify(arrayOfyVals))
			drawDiagram4Y( id, xVals, 'm', arrayOfyVals, yTitles, title, bottomLines, minMaxDefault);
		}
		this.readYVals = function(items, addDeltaT, jFrom) {
			var j, yVals = [];
			if (typeof jFrom === 'undefined') jFrom = 0;
			for (j = jFrom; j < items.length; j++) yVals.push(items[j].tRock + (addDeltaT && items[j].deltaT ? items[j].deltaT : 0));
//if (yVals.length == 0) throw new Error();
			return yVals;
		}
		this.drawIceStorage = function (ctx, dg) {
console.log("this.drawIceStorage...");
			const hIce = metrics.externalCryoStorage ? this.cycleData.iceStorage.height : 0;
			ctx.beginPath();
			//ctx.moveTo(toX(dg,0),toY(dg,this.tMinIceStorage)); 
			//ctx.moveTo(toX(dg,hIce),toY(dg,this.tMinIceStorage)); 
			ctx.moveTo(toX(dg,hIce),toY(dg,this.tMaxIceStorage)); 
			ctx.lineTo(toX(dg,0),toY(dg,this.tMaxIceStorage)); 
			ctx.lineTo(toX(dg,0),toY(dg,this.tMinIceStorage)); 
			//ctx.closePath();
			ctx.strokeStyle = 'blue';
			ctx.lineWidth = 2;
			ctx.stroke();
			ctx.font = "12px Arial";
			ctx.fillStyle = 'blue';
			var text = "tMax="+Math.round(this.tMaxIceStorage)+" K";
			ctx.fillText( text, toX(dg,hIce/2),toY(dg,this.tMaxIceStorage)); // -ctx.measureText(text).width
			ctx.font = "16px Arial";
			text = myRound(this.cycleData.iceStorage.totalMassInTons/1000000,1)+(metrics.externalIceStorage ? " mt crushed ice" : " mt crushed stone");
			ctx.fillText( text, toX(dg,0)+10,toY(dg,this.tMaxIceStorage/2));

			ctx.beginPath();
			//ctx.moveTo(toX(dg,hIce+10),toY(dg,this.tMinIceStorage-3));
			//ctx.moveTo(toX(dg,hLP),toY(dg,this.tMinIceStorage));
			ctx.moveTo(toX(dg,hLP),toY(dg,this.tMax)-3);
			ctx.lineTo(toX(dg,hIce),toY(dg,this.tMax)-3);
			ctx.lineTo(toX(dg,hIce),toY(dg,this.tMinIceStorage));
			//ctx.closePath();
			ctx.strokeStyle = 'black';
			ctx.lineWidth = 2;
			ctx.stroke();
			text = myRound(metrics.lowPressureStorageVolume * metrics.lpGravelDensity/1000000,1)+" Mt crushed stone";
			ctx.fillStyle = 'black';
			ctx.fillText( text, toX(dg,hIce)+60,toY(dg,this.tMax/2));

			ctx.font = "12px Arial";
			ctx.fillStyle = 'black';
			const iLast = lpArraysOfVals.length-1;
			const lpHeight = this.stepLen * this.lpTemps.length;
			text = "+"+myRound(lpArraysOfVals[iLast-1][0]-lpArraysOfVals[iLast][0],1)+" K";
			ctx.fillText( text, x=toX(dg,0) - ctx.measureText(text).width, y=toY(dg,lpArraysOfVals[iLast-1][0]));
//console.log("iLast="+iLast+", lpHeight="+lpHeight+", obj:"+JSON.stringify(lpArraysOfVals[iLast-1][0]));
//console.log("LP exits:"+x+","+y+":"+text);
			text = myRound(lpArraysOfVals[iLast][lpArraysOfVals[iLast].length-1]-lpArraysOfVals[iLast-1][lpArraysOfVals[iLast-1].length-1],0)+" K";
			ctx.fillText( text, x=toX(dg,lpHeight), y=toY(dg,lpArraysOfVals[iLast][lpArraysOfVals[iLast].length-1]));
//console.log("LP exits:"+x+","+y+":"+text);

		}
		this.enthalpyLookup = function( aH, t ) {
			const i = parseInt(Math.floor(t));
			const part = t - i;
			const h = aH[i]*(1 - part) + aH[i+1]*part;
			if (Number.isNaN(h)) throw new Error();
			return h;
		};
		this.getMaxWetnessData = function( wetGasMax, gasMaxWeight )
		{
			// gravel pressure control tank works best in a higher pressure but the max pressure of ice must be the low pressure
			const pTop = this.lpHeCCache.p * 1.35; // (metrics.externalIceStorage ? 1 : 1.35); // getFloatParam('pControlTopP', this.lpHeCCache.p/1000)*1000;
			const hLiq = Module.PropsSI('H', 'P', pTop, 'Q', 0, name);
			const hGas = Module.PropsSI('H', 'P', pTop, 'Q', 1, name);
			const tMax = Module.PropsSI('T', 'P', pTop, 'Q', 1, name);
			const t0 = Math.floor(this.tLowLp-2);
			const tMin = this.tLowLp;
			const hTP = Module.PropsSI('H', 'T', this.tLowLp, 'Q', 1, name);
			var rockH = [], iceH = []
			//const tMin = Module.PropsSI('T', 'H', hLiq, 'Q', 1, name);
			var h;
			//const dLiq = Module.PropsSI('D', 'P', this.lpHeCache.p, 'T', tMax-0.1, name);
			if (iceH.length == 0) {
				var t, h = 0;
				for (t = t0; t <= this.tHighHp+5; t++) {
					const cpRock = getRockHeatCapasityJinKgK(t);	// per l
					h = h + cpRock;		// per m3
					rockH.push(h);
				}
				if (metrics.externalIceStorage) {
					h = 0;
					for (t = t0; t <= 273; t++) {
						const cp = getIceHeatCapasityKJinKgK(t)*1000;	// ice density
						h = h + cp;	// per m3
						iceH.push(h);
					}
				}
				else {
					iceH = rockH; // no ice anywhere
				}
	//console.log(JSON.stringify(iceH));
			}
			this.enthalpyLookup0 = this.enthalpyLookup(metrics.externalIceStorage ? iceH : rockH, this.tLowLp-t0);
			const hDiff = this.enthalpyLookup(metrics.externalIceStorage ? iceH : rockH, tMax-t0) - this.enthalpyLookup0;
			const maxWetness = hDiff/(hGas-hLiq);
			this.pControlVolume = wetGasMax/(hDiff * (metrics.externalIceStorage ? metrics.crushedIceDensity : metrics.lpGravelDensity)*1000 / (hGas-hLiq));
			metrics.pControlVolume = this.pControlVolume;
			if (this.cycleData.iceStorage) {
				this.cycleData.iceStorage.pControlVolume = this.pControlVolume;
				metrics.rPControl = this.rPControl = Math.sqrt(this.pControlVolume/30/Math.PI);	// 20 m ice on top => max 250 kpa
				metrics.hPressureControl = 30;
				this.pControlArea = Math.PI*this.rPControl**2;
				this.cycleData.iceStorage.rPControl = this.rPControl;
				console.log("rPControl="+this.rPControl+", this.pControlArea="+this.pControlArea+", this.pControlVolume="+this.pControlVolume);
			}
			else if (metrics.pControlWithSolids) {
				// external pressure control storage of crushed stone  (pressurized to about 200 kPa)
				metrics.rPControl = this.rPControl = Math.sqrt(metrics.rGroundBelowActual**2/3);
				this.pControlArea = Math.PI*this.rPControl**2;
				metrics.hPressureControl = metrics.pControlVolume/this.pControlArea;
				if (metrics.hPressureControl < 25) {
					metrics.rPControl = this.rPControl = Math.sqrt(this.pControlVolume/25/Math.PI);	// 20 m ice on top => max 250 kpa
					metrics.hPressureControl = 25;
					this.pControlArea = Math.PI*this.rPControl**2;
				}
				console.log("rPControl="+this.rPControl+", this.pControlVolume="+this.pControlVolume+", metrics.hPressureControl="+metrics.hPressureControl);
			}
			//const hRockSat = this.enthalpyLookup(rockH, tMax-t0);
			var ret = {};
			//ret.height = this.hIce;
			ret.capasityKg = wetGasMax;
			//maxWetness *= 0.1;	
			ret.gasMaxWeight = gasMaxWeight;
			ret.maxWetness = maxWetness;
			ret.hMin = hTP;
			ret.hLiq = hLiq;
			ret.hGas = hGas;
			ret.pMax = pTop;
			ret.tMax = tMax;
			ret.t0 = t0;
			ret.tMin = tMin;
			ret.pMin = Module.PropsSI('P', 'T', this.tLowLp, 'Q', 1, name);
			ret.enthalpyLookup0 = this.enthalpyLookup0;
console.log("**** maxWetness: "+maxWetness+" => tMin="+tMin+", ret="+JSON.stringify(ret));
			ret.rockH = rockH;
			ret.iceH = iceH;
			return ret;
		};
		this.getWetnessSnapshot = function( title, maxWetnessData, pLow, pHigh, pHighCharge )
		{
			const t0 =  maxWetnessData.t0;
			const tMax =  maxWetnessData.tMax;
			var k, d, d2, h = this.stepLen/2
			var weightDiffKg = 0
			var tPrev = 0;
			for (k = 0; k < this.lpTemps.length; k++) {
				if (Math.abs(this.lpTemps[k].tRock - tPrev)>1)
				{
					d = Module.PropsSI('D', 'P', this.lpHeCache.p, 'T', this.lpTemps[k].tRock, name);
					d2 = this.lpHeCache.p == this.lpHeCCache.p ? d : Module.PropsSI('D', 'P', this.lpHeCCache.p, 'T', this.lpTemps[k].tRock, name);
					tPrev = this.lpTemps[k].tRock;
				}
				if (h < this.hIce) {
					weightDiffKg += (d-d2)*this.lpTemps[k].area*this.stepLen*this.crushedIceOrRockPorosity;
				}
				else {
					weightDiffKg += (d-d2)*this.lpTemps[k].area*this.stepLen*metrics.lpGravelPorosity;
				}
				h += this.stepLen;
			}
			var ret = {};
			ret.curChangeKg = Math.abs(weightDiffKg);
			tPrev = 0;
			for (k = 0; k < this.hpTemps.length; k++) {
				if (Math.abs(this.hpTemps[k].tRock - tPrev)>1)
				{
					d2 = Module.PropsSI('D', 'P', pHigh, 'T', this.hpTemps[k].tRock, name);
					d = Module.PropsSI('D', 'P', pHighCharge, 'T', this.hpTemps[k].tRock, name);
					tPrev = this.hpTemps[k].tRock;
				}
				weightDiffKg += (d2-d)*this.hpTemps[k].area*this.stepLen*metrics.hpGravelPorosity;
			}
			var tFrom, pFrom, hLiq, hGas, hDiff, wetGasKg;
			if (!testResult( weightDiffKg, t0, tMax, 50, function (t) {
				tFrom = t;
				hLiq = Module.PropsSI('H', 'T', t, 'Q', 0, name);
				hGas = Module.PropsSI('H', 'T', t, 'Q', 1, name);
				hDiff = this.enthalpyLookup(metrics.externalIceStorage ? maxWetnessData.iceH : maxWetnessData.rockH, t-t0) - maxWetnessData.enthalpyLookup0;
				curWetness = hDiff/(hGas-hLiq);
				wetGasKg = metrics.pControlVolume * (hDiff * (metrics.externalIceStorage ? metrics.crushedIceDensity : metrics.lpGravelDensity)*1000 / (hGas-hLiq));
				return wetGasKg;
				}, 20)) 
			{
console.log("t="+tFrom+", hDiff="+hDiff+", (hGas-hLiq)="+(hGas-hLiq)+" => wetGasKg="+wetGasKg+", weightDiffKg="+weightDiffKg);
throw new Error();
			}
			pFrom = Module.PropsSI('P', 'T', tFrom, 'Q', 1, name);

			ret.wetness = curWetness; //(hCurGas - maxWetnessData.h)/(hCurGas - hCurLiq);
			ret.t = tFrom;
			ret.p = pFrom;
			ret.title = title;
			ret.curKapasityKg = weightDiffKg;
console.log(title+" snapshot: "+JSON.stringify(ret));
			return ret;
		};
		this.fnSimulateChargeDischarge = function(noUI, simulations) {
			var lpArraysOfVals = [];
			var hpArraysOfVals = [];
			var bottomLines = [];
			var maxLpOutDeltaT = 0, maxHpOutDeltaT = 0;
			var minLpOutDeltaT = 0, minHpOutDeltaT = 0;
			var hTotal = 0;
			
			if (name == 'Hydrogen' || name == 'Helium') {
				alert("Storage simulatio is not possible with Hydrogen or Helium. They have too low saturation temperatures.");
				modal.style.display = "none";
				return;
			}
			var startTime = Date.now();
			//document.getElementById("modalStartButton").style.display = "none";
			//document.getElementById("modalStopButton").style.display = "block";
			var i, k, hLastDiscarge,  hCount;
			if (!simulations) simulations = getIntParam("heSimulations", 2 );

			for (i = 0; i < simulations; i++) {
				hCount = 0;
				var tLpOut, tHpOut;
				var prevCheck = -1;
				//modTitle.innerText = "Discharging the storage ...";
				//modContent.innerText = "hour 0";
				//sleep(50);

				if (noUI && this.gasWeightMaxDiff) {
					// gasWeight: empty storage in discharge mode, this.gasLpWeightMin: charged storage in charge mode
					maxWetnessData = this.getMaxWetnessData( metrics.hasLiquidStorage ? this.gasWeightMaxDiff : this.gasMaxWeight - this.gasLpWeightMin, this.gasMaxWeight );
				}
				do {
					if (!this.simulateDischarge( calcHeatExchangeCoeffcient, 1)) {
						console.log("hCount="+hCount);
						return;
					}
					//console.log(hCount+":9.t="+this.hpTemps[9].tRock+",h="+this.hpTemps[9].hGasOut+";10.t="+this.hpTemps[10].tRock+",h="+this.hpTemps[10].hGasOut+";11.t="+this.hpTemps[11].tRock+",h="+this.hpTemps[11].hGasOut);
					hCount += this.rockHeatingTime/3600;
					/*if (metrics.pControlWithSolids) {
						if (i == 1 && prevCheck != Math.round(hCount/50)) {
							prevCheck = Math.round(hCount/50);
							const ret = this.getWetnessSnapshot( "Discharge "+hCount, maxWetnessData, this.pLow, this.pHigh );
							ret.hour = hCount;
							disthis.chargeSnapshots.push(ret);
						}
					}*/
					if (Math.abs(minLpOutDeltaT) < Math.abs(this.lpTemps[this.hFrom].deltaT)) minLpOutDeltaT = this.lpTemps[this.hFrom].deltaT;
					if (Math.abs(maxHpOutDeltaT) < Math.abs(this.hpTemps[this.hpTemps.length-1].deltaT)) maxHpOutDeltaT = this.hpTemps[this.hpTemps.length-1].deltaT;
					tLpOut = this.lpTemps[this.hFrom].tRock + this.lpTemps[this.hFrom].deltaT;
					tHpOut = this.hpTemps[this.hpTemps.length-1].tRock + this.hpTemps[this.hpTemps.length-1].deltaT;
				} while (tHpOut>(this.tHighHp-this.hpMaxHighDiff) && tLpOut<(this.tLowLp+this.lpMaxLowDiff) && !this.stopChargeDischargeSimulation);
				/*if (i == 1 && metrics.pControlWithSolids)
				{
					const ret = this.getWetnessSnapshot( "Discharge "+hCount, maxWetnessData, this.pLow, this.pHigh );
					ret.hour = hCount;
					disthis.chargeSnapshots.push(ret);
				}*/

				if (typeof this.midHour === 'undefined') {
					this.midHour = parseInt(Math.round(hCount/2*this.cycleData.hotCycle[1].massFlow/this.cycleData.hotCycle[0].massFlow));
				}
				// calculate the weight of gas after discharge when all gas is in the storage
				for (k = 0; k < this.hpTemps.length; k++) {
					this.hpTemps[k].tMin = this.hpTemps[k].tRock;
				}
				if (!this.gasMaxWeight &&  !this.stopChargeDischargeSimulation)
				{
					// calc total gas weight in discarge
					this.gasHpWeight = 0;
					this.gasMaxWeight = 0, this.gasWeightMaxDiff = 0, this.gasLpWeight = 0, this.gasWeightMinDiff = 0; 
					var d, d2, tPrev = 0, h = this.stepLen/2
					for (k = 0; k < this.lpTemps.length; k++) {
						if (Math.abs(this.lpTemps[k].tRock - tPrev)>1)
						{
							d = Module.PropsSI('D', 'P', this.lpHeCache.p, 'T', this.lpTemps[k].tRock, name);
							d2 = this.lpHeCache.p == this.lpHeCCache.p ? d : Module.PropsSI('D', 'P', this.lpHeCCache.p, 'T', this.lpTemps[k].tRock, name);
							tPrev = this.lpTemps[k].tRock;
						}
						if (h < this.hIce) {
							this.gasLpWeight += d*this.lpTemps[k].area*this.stepLen*this.crushedIceOrRockPorosity;
							this.gasWeightMin2 += d2*this.lpTemps[k].area*this.stepLen*this.crushedIceOrRockPorosity;
							this.gasWeightMaxDiff += (d-d2)*this.lpTemps[k].area*this.stepLen*this.crushedIceOrRockPorosity;
						}
						else {
							this.gasLpWeight += d*this.lpTemps[k].area*this.stepLen*metrics.lpGravelPorosity;
							this.gasWeightMin2 += d2*this.lpTemps[k].area*this.stepLen*metrics.lpGravelPorosity;
							this.gasWeightMaxDiff += (d-d2)*this.lpTemps[k].area*this.stepLen*metrics.lpGravelPorosity;
						}
						h += this.stepLen;
					}
					if (this.gasWeightMaxDiff) {
						console.log("LP difference="+this.gasWeightMaxDiff+"???");
					}
					tPrev = 0;
					for (k = 0; k < this.hpTemps.length; k++) {
						if (Math.abs(this.hpTemps[k].tRock - tPrev)>1)
						{
							d = Module.PropsSI('D', 'P', this.hpHeCache.p, 'T', this.hpTemps[k].tRock, name);
							d2 = Module.PropsSI('D', 'P', this.hpHeCCache.p, 'T', this.hpTemps[k].tRock, name);
							tPrev = this.hpTemps[k].tRock;
						}
						this.gasHpWeight += d*this.hpTemps[k].area*this.stepLen*metrics.hpGravelPorosity;
						this.gasWeightMin2 += d2*this.hpTemps[k].area*this.stepLen*metrics.hpGravelPorosity;
					}
					this.gasMaxWeight = this.gasLpWeight + this.gasHpWeight;
					this.gasWeightMaxDiff = this.gasMaxWeight - this.gasWeightMin2;
console.log("gasWeightMaxDiff="+this.gasWeightMaxDiff+", this.gasMaxWeight="+this.gasMaxWeight);
				}
				hTotal += hCount;
				console.log("Reason NOT:"+tHpOut+">"+(this.tHighHp-this.hpMaxHighDiff)+" && "+tLpOut+"<"+(this.tLowLp+this.lpMaxLowDiff)+", hCount="+hCount);
	//console.log("hpLiquidStorage="+JSON.stringify(hpLiquidStorage));
				console.log("Discharged "+hCount+" hours, minLpOutDeltaT="+minLpOutDeltaT+", maxHpOutDeltaT="+maxHpOutDeltaT);
				lpArraysOfVals.push(this.readYVals(this.lpTemps, false, this.hFrom));
				hpArraysOfVals.push(this.readYVals(this.hpTemps));
				bottomLines.push((2*i+1)+". gravel temp after "+Math.round(hCount)+" hours discharge"); // : "+getUnusedLength(this.lpHeCache));
				hLastDiscarge = hCount;
				if (calcHeatExchangeCoeffcient && i == (simulations-1)) {
					lpArraysOfVals.push(this.readYVals(this.lpTemps, true, this.hFrom));
					hpArraysOfVals.push(this.readYVals(this.hpTemps, true));
					bottomLines.push("Gas temperature in last discharge (dots)"); // : "+getUnusedLength(this.lpHeCCache));
				}
//console.log("noUI="+noUI+", maxWetnessData="+JSON.stringify(maxWetnessData));
				//modTitle.innerText = "Charging the storage ...";
				//modContent.innerText = "hour 0";
				//sleep(50);
				hCount = 0;
				prevCheck = -1;
				for (;;) {
					if (!this.simulateCharge( calcHeatExchangeCoeffcient, 1)) return;
					if (false && i == 0 && hCount==0) {
						console.log(JSON.stringify(this.hpTemps));
						return;
					}
					hCount += this.rockHeatingTime/3600;

					if (i == 1 && prevCheck != Math.round(hCount/50)) {
						prevCheck = Math.round(hCount/50);
						const ret = this.getWetnessSnapshot( "Charge "+hCount, maxWetnessData, this.pLow, this.pHigh, this.hpHeCCache.p );
						ret.hour = hCount;
						this.chargeSnapshots.push(ret);
					}
					if (Math.abs(maxLpOutDeltaT) < Math.abs(this.lpTemps[this.lpTemps.length-1].deltaT)) maxLpOutDeltaT = this.lpTemps[this.lpTemps.length-1].deltaT;
					if (Math.abs(minHpOutDeltaT) < Math.abs(this.hpTemps[0].deltaT)) minHpOutDeltaT = this.hpTemps[0].deltaT;

					tLpOut = this.lpTemps[this.lpTemps.length-1].tRock + this.lpTemps[this.lpTemps.length-1].deltaT;
					tHpOut = this.hpTemps[0].tRock + this.hpTemps[0].deltaT;
					if (tHpOut>(this.tLowHp+this.hpMaxLowDiff) || tLpOut<=(this.tHighLp-this.lpMaxHighDiff) || this.stopChargeDischargeSimulation) 
					{
						console.log("Reason NOT:"+tLpOut+">"+(this.tHighLp-this.lpMaxHighDiff)+" && "+tHpOut+"<"+(this.tLowHp+this.hpMaxLowDiff)+", hCount="+hCount);
						break;
					}
					/*if (!swaptHighLpAndtLowHp) {
						if (tLpOut<=(this.tHighLp-this.lpMaxHighDiff)) {
							console.log("Reason:"+tLpOut+">"+(this.tLowHp+this.hpMaxLowDiff)+", hCount="+hCount);
							break;
						}
					}*/
				}
				for (k = 0; k < this.hpTemps.length; k++) {
					this.hpTemps[k].tMax = this.hpTemps[k].tRock;
				}
				if (i == 1) {
					const ret = this.getWetnessSnapshot( "Charge "+hCount, maxWetnessData, this.pLow, this.pHigh, this.hpHeCCache.p );
					ret.hour = hCount;
					this.chargeSnapshots.push(ret);
				}
				if (!this.gasWeightMin && !this.stopChargeDischargeSimulation)
				{
					// calc total gas weight in discarge
					this.gasWeightMin = 0, this.gasLpWeightMin = 0, this.gasHpWeightMin = 0;
					this.gasWeightMinDiff = 0; 
					var d, d2, tPrev = 0, h = this.stepLen/2
					for (k = 0; k < this.lpTemps.length; k++) {
						if (Math.abs(this.lpTemps[k].tRock - tPrev)>1)
						{
							d = Module.PropsSI('D', 'P', this.lpHeCache.p, 'T', this.lpTemps[k].tRock, name);
							d2 = Module.PropsSI('D', 'P', this.lpHeCCache.p, 'T', this.lpTemps[k].tRock, name);
							tPrev = this.lpTemps[k].tRock;
						}
						if (h < this.hIce) {
							this.lpGasWeightMinDiff += (d-d2)*this.lpTemps[k].area*this.stepLen*this.crushedIceOrRockPorosity;
							this.gasLpWeightMin += d2*this.lpTemps[k].area*this.stepLen*this.crushedIceOrRockPorosity;
							this.gasWeight2 = d*this.lpTemps[k].area*this.stepLen*this.crushedIceOrRockPorosity;
						}
						else {
							this.lpGasWeightMinDiff += (d-d2)*this.lpTemps[k].area*this.stepLen*metrics.lpGravelPorosity;
							this.gasLpWeightMin += d2*this.lpTemps[k].area*this.stepLen*metrics.lpGravelPorosity;
							this.gasWeight2 = d*this.lpTemps[k].area*this.stepLen*metrics.lpGravelPorosity;
						}
						h += this.stepLen;
					}
					h = this.stepLen/2, tPrev = 0;
					for (k = 0; k < this.hpTemps.length; k++) {
						if (Math.abs(this.hpTemps[k].tRock - tPrev)>1)
						{
							d = Module.PropsSI('D', 'P', this.hpHeCache.p, 'T', this.hpTemps[k].tRock, name);
							d2 = Module.PropsSI('D', 'P', this.hpHeCCache.p, 'T', this.hpTemps[k].tRock, name);
							tPrev = this.hpTemps[k].tRock;
						}
						this.gasWeightMinDiff += (d-d2)*this.hpTemps[k].area*this.stepLen*metrics.hpGravelPorosity;
						this.gasHpWeightMin += d2*this.hpTemps[k].area*this.stepLen*metrics.hpGravelPorosity;
						this.gasWeight2 = d*this.hpTemps[k].area*this.stepLen*metrics.hpGravelPorosity;
					}
					this.gasWeightMin = this.gasHpWeightMin + this.gasLpWeightMin;
console.log("gasLpWeightMin + this.gasHpWeightMin: "+this.gasLpWeightMin+" + "+this.gasHpWeightMin+"="+this.gasWeightMin+", this.gasWeightMinDiff="+this.gasWeightMinDiff+", this.hpHeCache.p="+this.hpHeCache.p+", this.hpHeCCache.p="+this.hpHeCCache.p);

					//this.getMaxWetnessData( gasWeight - this.gasWeightMin );

/*					
					var lpGasWeightMinDiffActual; 
console.log("discharge ratio: "+this.hpHeCache.p+"/"+this.lpHeCCache.pLowDischarge+"="+this.hpHeCache.p/this.lpHeCCache.pLowDischarge+", pLowCharge="+this.lpHeCCache.p+", pLowDischarge="+this.lpHeCCache.pLowDischarge);
					testResult( this.hpHeCache.p/this.lpHeCCache.pLowDischarge, this.lpHeCCache.p, this.lpHeCCache.pLowDischarge**2/this.lpHeCCache.p, 0.005, function(p) {
						pLowDischargeActual = p;
						this.lpGasWeightMinDiffActual = lpGasWeightMinDiff*(this.lpHeCCache.p - p) / (this.lpHeCCache.p - this.lpHeCCache.pLowDischarge);
						pHighDischargeActual = this.hpHeCCache.p * (this.gasHpWeightMin + lpGasWeightMinDiffActual) / this.gasHpWeightMin;
if (Number.isNaN(pLowDischargeActual) || Number.isNaN(lpGasWeightMinDiffActual) || Number.isNaN(pHighDischargeActual)) throw Error("Debug bp");
						return pHighDischargeActual/pLowDischargeActual;
					}, 15);
console.log("New discharge ratio: "+pHighDischargeActual+"/"+pLowDischargeActual+"="+pHighDischargeActual/pLowDischargeActual);
					this.lpGasWeightMinDiff = lpGasWeightMinDiffActual;
*/
					if (noUI || typeof maxWetnessData === 'undefined') {
						var maxDiff = this.gasWeightMaxDiff;
						var maxWeight = this.gasMaxWeight;
						var minWeight = this.gasWeightMin;
						if (this.gasWeight2 > this.gasMaxWeight) this.gasMaxWeight = this.gasWeight2;
						if (minWeight > this.gasWeightMin2) minWeight = this.gasWeightMin2;
						if (maxDiff < this.gasWeightMinDiff) maxDiff = this.gasWeightMinDiff;
						if (maxDiff < (maxWeight - minWeight)) maxDiff = (maxWeight - minWeight);
						if (metrics.hasLiquidStorage) maxDiff = this.gasWeightMaxDiff;
						metrics.maxWetnessData = maxWetnessData = 
							metrics.pControlWithSolids ? this.getMaxWetnessData( maxDiff, maxWeight ) : getNewWetnessData( maxDiff, maxWeight );
					}
					{
						var hpCapasityJ = 0;
						const tMin = maxWetnessData.t0;
						for (k = 0; k < this.hpTemps.length; k++) {
							const hMin = this.enthalpyLookup(maxWetnessData.rockH, this.hpTemps[k].tMin - maxWetnessData.t0) - maxWetnessData.enthalpyLookup0;
							const hMax = this.enthalpyLookup(maxWetnessData.rockH, this.hpTemps[k].tMax - maxWetnessData.t0) - maxWetnessData.enthalpyLookup0;
							hpCapasityJ += (hMax-hMin)*this.hpTemps[k].area*this.stepLen*metrics.gravelDensity*1000;
if (Number.isNaN(hpCapasityJ)) throw new Error();
						}
						this.cycleData.hpStorageHeatCapacity = hpCapasityJ;
						metrics.hpStorageUtilization = this.cycleData.hpStorageHeatCapacity/(1000*getAverageRockHeatCapasitykJinKg(this.cycleData.coldCycle[0].tOut,this.cycleData.hotCycle[0].tIn)*metrics.highPressureStorageVolume*1000*metrics.gravelDensity*(this.cycleData.hotCycle[0].tIn-this.cycleData.coldCycle[0].tOut));
						$("#hpStorageUtilization").val(myRound(metrics.hpStorageUtilization*100,1));
						console.log("hpCapasityJ=hpStorageHeatCapacity="+this.cycleData.hpStorageHeatCapacity+" J == "+Math.round(this.cycleData.hpStorageHeatCapacity/(3600*1000000))+" MWh, hpStorageUtilization="+metrics.hpStorageUtilization);
					}
					// densities OK, but otherwise these do not make very much sense (volumens OK)
console.log("HP: this.gasWeightMaxDiff tons="+this.gasWeightMaxDiff/1000+", this.gasWeightMinDiff tons="+this.gasWeightMinDiff/1000);
					const tGasMaxDischarge = Module.PropsSI('T', 'P', this.hpHeCache.p, 'Q', 1, name)+0.1;
					//const hGasMaxDischarge = Module.PropsSI('H', 'P', this.hpHeCache.p, 'T', tGasMaxDischarge, name);
					const dGasMaxCharge2 = Module.PropsSI('D', 'P', this.hpHeCCache.p, 'Q', 1, name)-0.1;
					const dGasMaxCharge = Module.PropsSI('D', 'P', this.hpHeCCache.p, 'T', tGasMaxDischarge, name);
console.log("max: dGasMaxCharge2="+dGasMaxCharge2+", possible: dGasMaxCharge="+dGasMaxCharge);
					const hGasMaxCharge = Module.PropsSI('H', 'P', this.hpHeCCache.p, 'D', dGasMaxCharge, name);
					const tGasMinCharge = Module.PropsSI('T', 'P', this.hpHeCCache.p, 'D', dGasMaxCharge, name);
					metrics.dCaveMin = Module.PropsSI('D', 'P', this.lpHeCCache.p, 'T', this.hpHeCCache.tHigh-100, name);	// guess the heat loss during the short fill period
					const dGasMaxDischarge = dGasMaxCharge*(this.gasMaxWeight-this.gasWeightMaxDiff)/this.gasMaxWeight;
					const hGasMaxDischarge = Module.PropsSI('H', 'P', this.hpHeCache.p, 'D', dGasMaxDischarge, name);
					const hGasColdHp = Module.PropsSI('H', 'P', this.hpHeCache.p, 'T', this.hpHeCCache.tLow, name);
					const mwhMax = (hGasMaxDischarge-hGasMaxCharge)*this.gasMaxWeight/1000000/3600;
					const hDeltaRock = getAverageRockHeatCapasitykJinKg(tGasMinCharge, this.hpHeCCache.tLow)*1000*(this.hpHeCCache.tLow-tGasMinCharge);
					const hDeltaGas = hGasColdHp - hGasMaxCharge;
					const coldRecoveryRockTons = this.gasWeightMaxDiff*hDeltaGas/hDeltaRock/1000;
					console.log("hDeltaRock="+hDeltaRock+", hDeltaGas="+hDeltaGas); 
					if (metrics.hPressureControlTank) {
						metrics.hCave = metrics.dome.hBottom*0.8+metrics.insCaveBottom;
						metrics.volCave = this.gasWeightMaxDiff/(dGasMaxCharge - metrics.dCaveMin);
						metrics.rCave = Math.sqrt(metrics.volCave/Math.PI/(metrics.hCave-metrics.insCaveBottom)) + metrics.insCaveBottom;
						metrics.dCaveMax = dGasMaxCharge;
						metrics.tCaveChargeMin = tGasMinCharge;
						metrics.pCaveDischargeMin = Module.PropsSI('P', 'D', metrics.dCaveMin, 'T', tGasMinCharge, name);
					}
if (Number.isNaN(metrics.hCave)) throw new Error();
console.log("metrics.hCave="+metrics.hCave+", metrics.dCaveMin="+metrics.dCaveMin+", tGasMinCharge="+tGasMinCharge);
console.log("dGasMax Charge/Discharge:"+dGasMaxCharge+"/"+dGasMaxDischarge+", mwhMax="+mwhMax+", coldRecoveryRockTons="+coldRecoveryRockTons+", case diameter="+metrics.rCave*2);
					const dGasMinCharge = dGasMaxCharge*this.gasWeightMin/this.gasMaxWeight;
					const hGasMinCharge = Module.PropsSI('H', 'P', this.hpHeCCache.p, 'D', dGasMinCharge, name);
					const dGasMinDischarge = dGasMinCharge*(this.gasWeightMin+this.gasWeightMinDiff)/this.gasWeightMin;
					const hGasMinDischarge = Module.PropsSI('H', 'P', this.hpHeCache.p, 'D', dGasMinDischarge, name);
					const mwhMin = (hGasMinDischarge-hGasMinCharge)*this.gasWeightMin/1000000/3600;
console.log("dGasMin Charge/Discharge:"+dGasMinCharge+"/"+dGasMinDischarge+", mwhMin="+mwhMin);

console.log("Min volume with liq storage: "+(this.gasWeightMaxDiff)/dGasMaxCharge );
console.log("Min volume without liq storage: "+(this.gasMaxWeight - this.gasWeightMin)/dGasMaxCharge);
console.log("Volume without liq storage + 50% pControl: "+(this.gasMaxWeight - this.gasWeightMin/2)/dGasMaxCharge);
					
					const dLiq = Module.PropsSI('D', 'P', this.lpHeCache.p, 'Q', 0, name);
					metrics.liqVolM3 = Math.abs(Math.round((this.gasMaxWeight-this.gasWeightMin)/dLiq));
					const gasWeightText = name+" min/max weight "+Math.round(this.gasWeightMin/1000)+"/"+Math.round(this.gasMaxWeight/1000)+" tons => max liquid volume: "+metrics.liqVolM3+" m3";
					var timeMin, timeMax;
					const switchTimesText = "Charge-to-discharge switch time with 100% power: from "+
						Math.round((timeMin = this.gasWeightMinDiff/1000/coldCycleG[0].massFlow)/60)+" to "+Math.round((timeMax = this.gasWeightMaxDiff/1000/coldCycleG[0].massFlow)/60) + " minutes";
if (!switchTimesText || typeof switchTimesText === 'undefined') throw new Error();
					const switchPowerMW = (coldCycleG[0].workConsumed+this.cycleData.heatPumpNetWorkIn-coldCycleG[1].workProduced)/1000;
					this.cycleData.switchMWhMin = switchPowerMW*timeMin*2/3600;
					this.cycleData.switchMWhMax = switchPowerMW*timeMax*2/3600;
console.log("switchWork MWh: "+this.cycleData.switchMWhMin+" - "+this.cycleData.switchMWhMax);
					/*
					const dischargePressureRange = "Discharge high/low pressure range (kPa):  "+Math.round(pHighDischargeActual/1000)+'/'+myRound(pLowDischargeActual/1000)+" - "
														+ Math.round(pHighDischargeActual2/1000)+'/'+myRound(pLowDischargeActual2/1000);
					console.log(dischargePressureRange);
					*/
					if (metrics.pressureControlInBottom && metrics.pHighBottomActual && typeof metrics.hPressureControlTank === 'undefined') {
						metrics.tPressureControlTank = this.cycleData.coldCycle[1].tIn;	// cooled to heat pump tLow
						const dMax = Module.PropsSI('D', 'T', metrics.tPressureControlTank, 'P', metrics.pHighBottomActual, name);
						const dMin = Module.PropsSI('D', 'T', metrics.tPressureControlTank, 'P', metrics.pLow, name);
						metrics.weigthPressureControlTank = metrics.hasLiquidStorage ? this.gasWeightMin : this.gasMaxWeight;
						const vol = (metrics.weigthPressureControlTank)/(dMax-dMin)/metrics.insGravelPorosity;	// max 0.42 
						metrics.volPressureControlTank = vol;
console.log("pressure control tank: dMax="+dMax+", 	gas weight="+metrics.weigthPressureControlTank+"=>"+", gas volume="+(metrics.weigthPressureControlTank)/(dMax-dMin));
						const rBottom = metrics.dome.rBottom;
						const maxAngleTan = Math.tan(Math.PI*metrics.dome.maxAngle/180);

						testResult( vol, rBottom, rBottom/2, vol*0.001, function (rPressureControlTank) {
							metrics.rPressureControlTank = rPressureControlTank;
							metrics.hPressureControlTank = (rBottom-rPressureControlTank)*maxAngleTan;
							return metrics.hPressureControlTank*Math.PI/3*(rBottom**2 + rBottom*rPressureControlTank + rPressureControlTank**2);
						}, 15);
						metrics.rPressureControlTankI = metrics.rPressureControlTank * metrics.dome.rTopI / metrics.dome.rTop;	// better insulation in the top?
						console.log("Storage max density="+dMax+", volume:"+vol+", height="+metrics.hPressureControlTank+", metrics.rPressureControlTank="+metrics.rPressureControlTank);
					}
					else if (typeof metrics.hPressureControlTank === 'undefined') {
						metrics.hPressureControlTank = 0;
					}
					if (metrics.hPressureControlTank) {
						metrics.hCave = 0;
						metrics.rCave = 0;
					}
					else if (false && !metrics.pControlWithSolids) {
						const volColdRecoveryRock = coldRecoveryRockTons/metrics.gravelDensity;
						// assume 45 degrees walls, h = rBottom, rTop = 2 * rBottom =>
						// vol = PI*(2*rBottom)**2*(2*rBottom)/3 - PI*rBottom*rBottom*rBottom/3
						// 	 = PI*2**3*rBottom**3/3 - PI*rBottom**3/3 = PI*7*rBottom**3/3
						// => rBottom = pow(3*vol/(PI*7), 1/3);
						metrics.rBottomHeatRecovery = Math.pow(3*volColdRecoveryRock/(Math.PI*7), 1/3);
						metrics.rBottomHeatRecoveryI = metrics.rBottomHeatRecovery - metrics.insCaveBottom;
						console.log("metrics.rBottomHeatRecovery="+metrics.rBottomHeatRecovery);
					}
					//document.getElementById("modGasWeight").innerText = gasWeightText;
					console.log(gasWeightText);
					console.log(switchTimesText);
					//console.log(dischargePressureRange);
					this.cycleData.totalGasWeightText = gasWeightText;
					this.cycleData.switchTimeText = switchTimesText;
if (typeof switchTimesText === 'undefined') throw new Error();
					if (noUI) return;
				}
			
				hTotal += hCount;
				console.log("Charged "+hCount+" hours, maxLpOutDeltaT="+maxLpOutDeltaT+", minHpOutDeltaT="+minHpOutDeltaT);
				lpArraysOfVals.push(this.readYVals(this.lpTemps, false, this.hFrom));
				hpArraysOfVals.push(this.readYVals(this.hpTemps));
				bottomLines.push((2*i+2)+". gravel temp after "+Math.round(hCount)+" hours charge"); // : "+getUnusedLength(this.lpHeCCache));
				
				if (calcHeatExchangeCoeffcient && i == (simulations-1)) {
					lpArraysOfVals.push(this.readYVals(this.lpTemps, true, this.hFrom));
					hpArraysOfVals.push(this.readYVals(this.hpTemps, true));
					bottomLines.push("Gas temperature in last charge (dots)"); // : "+getUnusedLength(this.lpHeCCache));
				}
				if (this.stopChargeDischargeSimulation) break;
			}
		    var elapsedTime = Date.now() - startTime;
			console.log("elapsedTime="+Math.round(elapsedTime/1000)+" seconds "); 
			
	console.log("Showing diagarm...");

			minMaxDefault = {min : Math.round(this.tLowLp/10)*10, max : (1+Math.round(this.tHighLp/10))*10};
			this.drawHeatExchangeDiagram2( 'lpHeatExchangeDiagram', this.lpTemps, Math.round(hTotal)+" hours * "+workToMW(this.cycleData.netWorkOut)+" MW simulation of low-pressure storage", 
				bottomLines, lpArraysOfVals, minMaxDefault, metrics.externalCryoStorage ? this.drawIceStorage : function (ctx, dg) {
									text = myRound(metrics.lowPressureStorageVolume * metrics.lpGravelDensity/1000000,1)+" Mt crushed stone";
									ctx.fillStyle = 'black';
									ctx.font = "16px Arial";
									ctx.fillText( text, toX(dg,height/4)+60,toY(dg,this.tHighLp/2));
									ctx.font = "12px Arial";
									const iLast = lpArraysOfVals.length-1;
									const lpHeight = this.stepLen * (this.lpTemps.length-this.hFrom);
									text = "+"+myRound(lpArraysOfVals[iLast-1][0]-lpArraysOfVals[iLast][0],1)+" K";
									ctx.fillText( text, x=toX(dg,0) - ctx.measureText(text).width, y=toY(dg,lpArraysOfVals[iLast-1][0]));
						//console.log("iLast="+iLast+", lpHeight="+lpHeight+", obj:"+JSON.stringify(lpArraysOfVals[iLast-1][0]));
						//console.log("LP exits:"+x+","+y+":"+text);
									text = Math.round(lpArraysOfVals[iLast][lpArraysOfVals[iLast].length-1]-lpArraysOfVals[iLast-1][lpArraysOfVals[iLast-1].length-1])+" K";
									ctx.fillText( text, x=toX(dg,lpHeight), y=toY(dg,lpArraysOfVals[iLast][lpArraysOfVals[iLast].length-1]));
						//console.log("LP exits:"+x+","+y+":"+text);
								});
			//var rect = document.getElementById('lpHeatExchangeDiagram').getBoundingClientRect();
			//console.log("lpHeatExchangeDiagram pos:"+rect.top, rect.right, rect.bottom, rect.left);
			this.drawStoneAndIceEnthalpies( 'usedHeatCapasitiesAndEnthalpies', this.tLowLp, this.tHighHp );
			minMaxDefault = {min : Math.round(this.tLowHp/100)*100, max : (1+Math.round(this.tHighLp/100))*100};
			this.drawHeatExchangeDiagram2('hpHeatExchangeDiagram', this.hpTemps, Math.round(hTotal)+ " hours * "+workToMW(this.cycleData.netWorkOut)+" simulation of high-pressure storage", bottomLines, hpArraysOfVals, minMaxDefault,
						function (ctx, dg) {
							text = myRound(metrics.highPressureStorageVolume * metrics.gravelDensity/1000000,1)+" mt crushed stone";
							//text2 = 'Storage capasity '+Math.round(hLastDiscarge*)+' GWh after '+Math.round(hTotal)+ " hours";
							ctx.font = "16px Arial";
							ctx.fillStyle = 'black';
							const height = metrics.dome.hBottomI;
//console.log("Hp storage text="+text+", pxX="+toX(dg,height/4)+", pxY="+toY(dg,this.tHighHp/2));							
							ctx.fillText( text, toX(dg,height/4),toY(dg,this.tHighHp/2));
							//ctx.fillText( text2, toX(dg,height/4),toY(dg,this.tHighHp/2)+20);
							ctx.font = "12px Arial";
							const iLast = hpArraysOfVals.length-1;
							text = "+"+myRound(hpArraysOfVals[iLast][0]-hpArraysOfVals[iLast-1][0],1)+" K";
							ctx.fillText(text, x=toX(dg,0)-ctx.measureText(text).width, y=toY(dg,hpArraysOfVals[iLast][0]));
//console.log("iLast="+iLast+", obj:"+JSON.stringify(hpArraysOfVals[iLast-1][0]));
//console.log("HP exits:"+x+","+y+":"+text);
							text = myRound(hpArraysOfVals[iLast-1][hpArraysOfVals[iLast-1].length-1] - hpArraysOfVals[iLast][hpArraysOfVals[iLast].length-1])+" K";
							ctx.fillText(text, x=toX(dg,height), y=toY(dg,hpArraysOfVals[iLast-1][hpArraysOfVals[iLast-1].length-1]));
//console.log("HP exits:"+x+","+y+":"+text+", t2="+hpArraysOfVals[iLast][hpArraysOfVals[iLast].length-1]+", t1="+hpArraysOfVals[iLast-1][hpArraysOfVals[iLast-1].length-1]);
						});
			this.drawPressureControlDiagram();
			this.drawRadiativeHeatTransferDiagram();
			
			const modal = document.getElementById("myModal");
			modal.style.display = "none";
			document.getElementById("lpHeatExchangeDiagram").scrollIntoView();
		}
		this.drawPressureControlDiagram = function() {
			var yVals = []; 
			var xVals = []; 
			var yTitles = [];
			var bottomLines = [];
			yVals.push([]); yVals.push([]); yVals.push([]); //yVals.push([]);yVals.push([]); yVals.push([]);
			xVals.push([]); xVals.push([]); xVals.push([]); //xVals.push([]);xVals.push([]); xVals.push([]);
			bottomLines.push("Tank pressure during charge");
			bottomLines.push("Tank temperature during charge");
			bottomLines.push(name+" wetness % during charge");
			yTitles.push("kPa");
			yTitles.push("K");
			yTitles.push("%");
//console.log("chargeSnapshots="+JSON.stringify(this.chargeSnapshots));
			var dischargeMid;
			var chargeMid;
			var lastBest = 100;
			var maxHour = this.chargeSnapshots[this.chargeSnapshots.length-1].hour;
			for (i = 0; i < this.chargeSnapshots.length; i++) {
				yVals[0].push( this.chargeSnapshots[i].p/1000 );
				yVals[1].push( this.chargeSnapshots[i].t );
				yVals[2].push( this.chargeSnapshots[i].wetness*100 );
				const percent = this.chargeSnapshots[i].hour / maxHour * 100;
				if (Math.abs(percent-50) < lastBest) {
					chargeMid = this.chargeSnapshots[i];
					lastBest = Math.abs(percent-50);
				}
				xVals[0].push(percent);
				xVals[1].push(percent);
				xVals[2].push(percent);
			}
			if (metrics.pControlVolume) {
				//console.log("Max/mid ratio="+(chargeMax.p/dischargeMid.p)+", pDischargeMid="+dischargeMid.p+" ("+(dischargeMid.hour / maxHour * 100)+"%), pChargeMid="+chargeMid.p+" ("+(chargeMid.hour / maxHour * 100)+"%)");
				minMaxDefault = {xMinMax:{min:0,max:100, step:10}, yMinMax:null};
				drawDiagram4Y( "lpStorageWetness", xVals, '%', yVals, yTitles, 
					"Pressure control tank temperature, pressure and wetness variation with storage charge level (during charge)", bottomLines, minMaxDefault,
					function (ctx, dg) {
								const text = metrics.externalIceStorage ?  
											myRound(metrics.pControlVolume * metrics.crushedIceDensity/1000000,1)+" Mt crushed ice" :
											myRound(metrics.pControlVolume * metrics.lpGravelDensity/1000000,1)+" Mt crushed stone";
									const text2 = "Radius: "+Math.round(metrics.rPControl)+" m, height: "+myRound(metrics.hPressureControl)+" m";
									const text3 = "Contains "+Math.round(this.chargeSnapshots[0].curKapasityKg/1000)+' - '+
													Math.round(this.chargeSnapshots[this.chargeSnapshots.length-1].curKapasityKg/1000)+' tons liquid '+name+' in charge';
													
									ctx.font = "16px Arial";
									ctx.fillStyle = 'black';
									ctx.fillText( text, dg.width/4, dg.height/2);
									ctx.fillText( text2, dg.width/4, dg.height/2+18);
									ctx.fillText( text3, dg.width/8, dg.height/2+36);
								});
			}
		}
		// radiative heat transfers in storage
		this.drawRadiativeHeatTransferDiagram = function()  {
			const mmMdpFrom = 0, mmMdpTo = 40, mmMdpStep = 1, areaFactor = 2, tDiffPerM = 200;
			var mdp, xVals = [], yVals = [], yTitles = [], bottomLines = [];
			for (mdp = mmMdpFrom; mdp <= mmMdpTo; mdp += mmMdpStep) xVals.push(mdp);
			yTitles.push('W/K/m');
			yVals.push(getRadiativeHeatTransfers( metrics.tHigh, mmMdpFrom, mmMdpTo, mmMdpStep, areaFactor, tDiffPerM ));
			bottomLines.push(Math.round(metrics.tHigh)+' K');
			yVals.push(getRadiativeHeatTransfers( metrics.tHigh-100, mmMdpFrom, mmMdpTo, mmMdpStep, areaFactor, tDiffPerM ));
			bottomLines.push(Math.round(metrics.tHigh-100)+' K');
			yVals.push(getRadiativeHeatTransfers( metrics.tHigh-200, mmMdpFrom, mmMdpTo, mmMdpStep, areaFactor, tDiffPerM ));
			bottomLines.push(Math.round(metrics.tHigh-200)+' K');
			yVals.push(getRadiativeHeatTransfers( metrics.tHigh-300, mmMdpFrom, mmMdpTo, mmMdpStep, areaFactor, tDiffPerM ));
			bottomLines.push(Math.round(metrics.tHigh-300)+' K');
			yVals.push(getRadiativeHeatTransfers( metrics.tHigh-400, mmMdpFrom, mmMdpTo, mmMdpStep, areaFactor, tDiffPerM ));
			bottomLines.push(Math.round(metrics.tHigh-400)+' K');
			yVals.push(getRadiativeHeatTransfers( this.tHighLp, mmMdpFrom, mmMdpTo, mmMdpStep, areaFactor, tDiffPerM ));
			bottomLines.push(Math.round(this.tHighLp)+' K');
			minMaxDefault = {xMinMax:{min:0,max:40, step:5}, yMinMax:null};

			drawRhtLine = function(ctx, dg, t, mmMdp, x, title) {
				const rht = getRadiativeHeatTransfer( t, mmMdp, areaFactor, tDiffPerM );
				const y = toY(dg, rht);
				const text = title +' ('+myRound(mmMdp)+' mm, '+Math.round(t)+' K) '+myRound(rht,2)+' W/K/m';
console.log(text+", y="+y);
				ctx.beginPath();
				ctx.moveTo(toX(dg,0), y);
				ctx.lineTo(x, y);
				ctx.stroke();
				ctx.fillText(text, x, y);
				return ctx.measureText(text).width;
			}
			drawDiagram4Y( "radiativeHeatTransfer", xVals, 'mm', yVals, yTitles, 
					"Radiative heat transfer of crushed stone and inslutations, variation with temperature and mean diameter of particles, temperature change 200 K/m; area factor 2", bottomLines, minMaxDefault,
					function (ctx, dg) {
						ctx.font = "12px Arial";
						ctx.fillStyle = 'black';
						ctx.strokeStyle = 'black';
						var mmMdp = 3/1000/0.07/0.42; // Stone wool fiber min diameter 3 um, min porosity 93% => 3/0.07 ¨ 0.043 mm
						var title = "Stone wool+screened crushed stone";
console.log(title+", mdp="+mmMdp);
						drawRhtLine(ctx, dg, this.hpHeCache.tHigh, mmMdp, toX(dg,5), title); 
						title = "Low-pressure storage";
						drawRhtLine(ctx, dg, this.lpHeCache.tHigh, this.lpHeCache.medianDP*1000, toX(dg,15), title); 
						title = "High-pressure storage";
						drawRhtLine(ctx, dg, this.hpHeCache.tHigh, this.hpHeCache.medianDP*1000, toX(dg,15), title); 
					});
		}
	}
	// end of constructor!

	firstChargeDischargeSimulation()
	{
		this.stopChargeDischargeSimulation = false;
		const noUI = true;
		this.fnSimulateChargeDischarge(noUI,2);
	}
	
	drawHeatExchangeDiagrams()
	{
		//const this.lpMaxLowDiff = 1, this.lpMaxHighDiff  = 20 , this.hpMaxLowDiff = 10, this.hpMaxHighDiff = 50;
console.log("Starting simulation...");
		this.stopChargeDischargeSimulation = false;
		const modal = document.getElementById("myModal");
		var lpGasWeightMinDiff=0;
		// the same pressure ratio & charge high pressure is the upper limit => 
		//var pLowInDischargeTEST = this.lpHeCache.p/this.lpHeCCache.p < this.hpHeCache.p/this.hpHeCCache.p ? this.hpHeCCache.p*this.lpHeCache.p/this.hpHeCache.p : this.lpHeCache.p;
		//const iceH = [], rockH = [];
		const modTitle = document.getElementById("modalTitle");
		const modContent = document.getElementById("modalContent");
		modTitle.innerText = "Press button to start the simulation";
		$(modContent).html("The charge/discharge simulation will take "+(simulations*30)+" seconds or more. The browser will be unresponsive. <bold>Press Control-Shift-I to see the progress in browser console and wait ...</bold>");
		modal.style.display = "block";
		document.getElementById("modalStopButton").style.display = "none";
		document.getElementById("modalStartButton").style.display = "block";
	}
}