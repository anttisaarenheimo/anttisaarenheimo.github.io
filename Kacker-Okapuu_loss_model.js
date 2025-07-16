// See: https://www.lth.se/fileadmin/tpe/Kurser/Chapter_4.pdf
function loss_model_KO(type, Re, angle_in, angle_out, c, s, H, t_cl, t_max, Ma_rel_in, Ma_rel_out, r_ht_in, p_in, p0rel_in, p_out, p0rel_out ) {
    // Compute the loss coefficient using the Kacker-Okapuu loss model
    // Author: Roberto Agromayor
	const t_te = s / 40;  // "For the T106C cascade, increasing the trailing edge thickness from 1.9% pitch to 2.8% pitch has a small effect on the loss"
	// incidence and deviation angles are assumed to be zero => metal angles are the same as angle_in and angle_out
	const staggerAngle = (angle_in + angle_out) / 2;	
	const b = c * Math.cos(staggerAngle);
	const o = s * Math.cos(angle_out);

    // Compute the loss coefficients
    // Mach number correction factor
    const f_Ma = mach_correction(Ma_rel_out);
//console.log(f_Ma+" = mach_correction("+Ma_rel_out+")");
    // const f_Ma = 1; // Uncomment this line to ignore the supersonic penalty

    // Reynolds number correction factor
    const f_Re = reynolds_correction(Re);
    // const f_Re = 1;
    // Profile loss coefficient
    let Y_p = profile_loss(type, angle_in, angle_out, c, s, t_max, 
        Ma_rel_in, Ma_rel_out, r_ht_in, p_in, p0rel_in, p_out, p0rel_out);

    // Corrected profile loss coefficient
    Y_p *= f_Ma * f_Re;

    // Secondary loss coefficient
    const Y_s = secondary_loss(angle_in, angle_out, Ma_rel_in, Ma_rel_out, H, c, b);

    // Clearance loss coefficient
    const Y_cl = clearance_loss(type, angle_in, angle_out, H, c, t_cl);

    // Trailing edge loss coefficient
    const Y_te = trailing_loss(angle_in, angle_out, t_te, o);

    // Overall loss coefficient
    const Y = Y_p + Y_s + Y_cl + Y_te;
    return { 'Y' : Y, 'Y_p' : Y_p, 'Y_s' : Y_s, 'Y_cl': Y_cl, 'Y_te' : Y_te };
}
function profile_loss(type, angle_in, angle_out, c, s, t_max, Ma_rel_in, Ma_rel_out, r_ht_in, p_in, p0rel_in, p_out, p0rel_out) {
    // Inlet shock loss
    let a = Math.max(0, f_hub(r_ht_in, type) * Ma_rel_in - 0.40);
//console.log(f_hub(r_ht_in, type)+" * "+Ma_rel_in+" - "+0.40+"="+a);
    let Y_shock = 0.75 * Math.pow(a, 1.75) * r_ht_in * (p0rel_in - p_in) / (p0rel_out - p_out);
//console.log("0.75 * "+Math.pow(a, 1.75)+" * "+r_ht_in+" * ("+p0rel_in+" - "+p_in+") / ("+p0rel_out+" - "+p_out+")="+Y_shock);
    // Avoid unphysical results if r_ht_in becomes negative during the optimization iterations
    Y_shock = Math.max(0, Y_shock);

    // Compressible flow correction factor
    // Limit excessively low values (it might be a problem during optimization)
    // The limit to 0.1 is quite arbitrary, but it worked for me
    let Kp = Math.max(0.1, K_p(Ma_rel_in, Ma_rel_out));

    // Compute the profile loss coefficient
    // The absolute value of the outlet angle is used as input
    // This is meaningless for stator cascades as angle_out > 0
    // However, for rotor cascades angle_out < 0 and the absolute value is used to
    // compute Y_p_reaction and Y_p_impulse according to Aungier correlation
    // These formulas are valid for 40 < abs(angle_out) < 80
    // Extrapolating outside of this limits might give completely wrong results
    // If the optimization algorithm has upper and lower bounds for the outlet 
    // angle there is no need to worry about this problem
    // angle_out_bis keeps the 40deg-losses for outlet angles lower than 40deg
    let angle_out_bis = Math.max(Math.abs(angle_out), 40 * Math.PI / 180);
    let Y_p_reaction = nozzle_blades(s / c, angle_out_bis);
    let Y_p_impulse = impulse_blades(s / c, angle_out_bis);

    // This formula works for both stator and rotor cascades
    // The logic for this formula is quite tricky and I had to put a lot of
    // thought into it. See my handwritten notes for information
    // There is no need to change the sign of any angle in this formula
    let Y_p = Y_p_reaction - Math.abs(angle_in / angle_out) * (angle_in / angle_out) * (Y_p_impulse - Y_p_reaction);

    // Limit the extrapolation of the profile loss to avoid negative values for
    // blade profiles with little deflection
    // Low limit to 80% of the axial entry nozzle profile loss
    // This value is completely arbitrary
    Y_p = Math.max(Y_p, 0.80 * Y_p_reaction);

    // Avoid unphysical effect on the thickness by defining the variable aa
    let aa = Math.max(0, -angle_in / angle_out);
    Y_p = Y_p * Math.pow((t_max / c) / 0.20, aa);
    Y_p = 0.914 * (2 / 3 * Y_p * Kp + Y_shock);
console.log("Y_p="+Y_p+" = 0.914 * (2 / 3 * "+Y_p+" * "+Kp+" + "+Y_shock+", Y_p_reaction="+Y_p_reaction+", Y_p_impulse="+Y_p_impulse+
		", Ma_rel_in="+Ma_rel_in+", Ma_rel_out="+Ma_rel_out+", angle_in/angle_out="+angle_in+"/"+angle_out+", (t_max / c)="+(t_max / c));
    return Y_p;
}
function secondary_loss(angle_in,angle_out,Ma_rel_in,Ma_rel_out,H,c,b)
{
	// Compute the secondary loss coefficient
	// Limit excessively low values (it might be a problem during optimization)
	// The limit to 0.1 is quite arbitrary, but it worked for me
	const Ks      = Math.max(0.10,K_s(Ma_rel_in,Ma_rel_out,H/b));
	const angle_m = Math.atan((Math.tan(angle_in)+Math.tan(angle_out))/2);
	const Z       = 4*(Math.tan(angle_in)-Math.tan(angle_out))**2*Math.cos(angle_out)**2/Math.cos(angle_m);
	const Y_s     = 1.2*Ks*0.0334*f_AR(H/c)*Math.cos(angle_out)/Math.cos(angle_in)*Z;
	return Y_s;
}


function clearance_loss(type, angle_in, angle_out, H, c, t_cl) {
    let B;
    if (type === 'stator') {
        B = 0.00; // Empirical parameter for the stator
		return 0;
    } else if (type === 'rotor') {
        B = 0.37; // Empirical parameter for the rotor (choose 0.37 for shrouded blades)
    } else {
        throw new Error("Specify the type of cascade: 'rotor' or 'stator'");
    }

    const angle_m = Math.atan((Math.tan(angle_in) + Math.tan(angle_out)) / 2);
    const Z = 4 * Math.pow((Math.tan(angle_in) - Math.tan(angle_out)), 2) * Math.pow(Math.cos(angle_out), 2) / Math.cos(angle_m);
    const Y_cl = B * (c / H) * Math.pow((t_cl / H), 0.78) * Z;

    return Y_cl;
}

function trailing_loss(angle_in, angle_out, t_te, o) {
    // Reacting blading
    const r_to_data_reaction = [0.000, 0.100, 0.200, 0.300, 0.400];
    const phi_data_reaction = [0.000, 0.017, 0.045, 0.090, 0.150];

    // Impulse blading
    const r_to_data_impulse = [0.000, 0.100, 0.200, 0.300, 0.400];
    const phi_data_impulse = [0.000, 0.010, 0.025, 0.047, 0.075];

    // Numerical trick to avoid too big r_to's
    const r_to = Math.min(0.400, t_te / o);

    // Interpolate data
    const d_phi2_reaction = linearInterpolation(r_to_data_reaction, phi_data_reaction, r_to);
    const d_phi2_impulse = linearInterpolation(r_to_data_impulse, phi_data_impulse, r_to);

    let d_phi2 = d_phi2_reaction - Math.abs(angle_in / angle_out) * (angle_in / angle_out) * (d_phi2_impulse - d_phi2_reaction);

    // Limit the extrapolation of the trailing edge loss
    d_phi2 = Math.max(d_phi2, d_phi2_impulse / 2);
    const Y_te = 1 / (1 - d_phi2) - 1;

    return Y_te;
}

function f_hub(r_ht, type) {
	/*
	Mach number of tip does not stop growing below hub/tip ratio 0.5!!
    if (r_ht < 0.50) {
        r_ht = 0.50; // Numerical trick to prevent extrapolation
    }
	*/

    // Stator curve
    const r_ht_data_S = [0.50, 0.60, 0.70, 0.80, 0.90, 1.00];
    const f_data_S = [1.40, 1.18, 1.05, 1.00, 1.00, 1.00];

    // Rotor curve
    const r_ht_data_R = [0.50, 0.60, 0.70, 0.80, 0.90, 1.00];
    const f_data_R = [2.15, 1.70, 1.35, 1.12, 1.00, 1.00];

    let f;
    if (type === 'stator') {
        f = linearInterpolation(r_ht_data_S, f_data_S, r_ht);
    } else if (type === 'rotor') {
        f = linearInterpolation(r_ht_data_R, f_data_R, r_ht);
    } else {
        throw new Error("Specify the type of cascade: 'rotor' or 'stator'");
    }

    return f;
}

function f_AR(x) {
    return (x < 2)  ? (1 - 0.25 * Math.sqrt(2 - x)) / x : 1 / x; //  (x >= 2);
}

function mach_correction(Ma_rel_out) {
    return 1 + ((Ma_rel_out > 1) ? Math.pow((Ma_rel_out - 1), 2)*60 : 0 );
}

function reynolds_correction(Re) {
	// Antti: why the correlation is exactly constant 1 between 2e5 and 1e6 ???
    return (Re < 2e5) ? Math.pow((Re / 2e5), -0.40) : (Re <= 1e6) ? 1 : Math.pow((Re / 1e6), -0.20); // (Re > 1e6);
}

function nozzle_blades(r_sc, angle_out) {
    // Use Aungier correlation to compute the pressure loss coefficient
    // This correlation is a formula that reproduces the figures from the Ainley
    // and Mathieson original figures

    const phi = 90 - angle_out * 180 / Math.PI;
    const r_sc_min = (phi < 30) ? (0.46 + phi / 77) : (0.614 + phi / 130); // (phi >= 30);
    const X = r_sc - r_sc_min;
	// Antti: The result is non continuous between out angle 27 and 30!!! Should 27 be 30???
    const A = (phi < 27) ? (0.025 + (27 - phi) / 530) : (0.025 + (27 - phi) / 3085);
    const B = 0.1583 - phi / 1640;
    const C = 0.08 * ((phi / 30) ** 2 - 1);
    const n = 1 + phi / 30;
    const Y_p_reaction = (phi < 30) ? (A + B * (X ** 2) + C * (X ** 3)) : (A + B * (Math.abs(X) ** n)); // (phi >= 30);

    return Y_p_reaction;
}

function impulse_blades(r_sc, angle_out) {
    // Use Aungier correlation to compute the pressure loss coefficient
    // This correlation is a formula that reproduces the figures from the Ainley
    // and Mathieson original figures

    const phi = 90 - angle_out * 180 / Math.PI;
    const r_sc_min = 0.224 + 1.575 * (phi / 90) - (phi / 90) ** 2;
    const X = r_sc - r_sc_min;
    const A = 0.242 - phi / 151 + (phi / 127) ** 2;
    const B = (phi < 30) ? (0.30 + (30 - phi) / 50)  : (0.30 + (30 - phi) / 275);
    const C = 0.88 - phi / 42.4 + (phi / 72.8) ** 2;
    const Y_p_impulse = A + B * (X ** 2) - C * (X ** 3);

    return Y_p_impulse;
}

function K_p(Ma_in, Ma_out) {
    return 1 - K_2(Ma_in / Ma_out) * (1 - K_1(Ma_out));
    // Kp = 1;
}

function K_s(Ma_in, Ma_out, r_Hb) {
    return 1 - K_3(1 / r_Hb) * (1 - K_p(Ma_in, Ma_out));
    // Ks = 1;
}

function K_1(x) {
    return (x < 0.20) ? 1 :  (x > 0.20 && x < 1.00)  ? + (1 - 1.25 * (x - 0.20)) : 0;
}

function K_2(x) {
    return x ** 2;
}

function K_3(x) {
    return x ** 2;
}

function linearInterpolation( xVals, yVals, x )
{
	if (x < xVals[0]) return yVals[0] + (yVals[1] - yVals[0])* (1 - (xVals[1]-x)/(xVals[1]-xVals[0]));
	var i;
	for (i = 0; i < xVals.length; i++) {
		if (x <= xVals[i+1] || (i+2) == xVals.length) {
			return yVals[i] + (yVals[i+1] - yVals[i])* (1 - (xVals[i+1]-x)/(xVals[i+1]-xVals[i]));
		}
	}
	throw new Error("getArrayCoefficient: i="+i+", xyVals="+JSON.stringify(xyVals)+", xVal="+xVal );
}

/*
function linearInterpolation(x, y, x0) {
    const idx = x.findIndex((val) => val >= x0);
    if (idx === -1 || idx === 0) {
        return y[0]; // extrapolate
    }
    const x1 = x[idx - 1];
    const x2 = x[idx];
    const y1 = y[idx - 1];
    const y2 = y[idx];
    return y1 + (y2 - y1) * ((x0 - x1) / (x2 - x1));
}
*/
