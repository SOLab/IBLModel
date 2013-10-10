import math
import cmath
from matplotlib import pyplot as plt


def lambdaCalc(epsilon, stratPar):
    ''' function calculate dimensionless function lambda
        epsilon - main fitting parameter of the model
        stratPar - stratification parameter, mu
    '''

    if stratPar > 0:
        lamb = (math.sqrt(1 + 40 * epsilon ** 2 * stratPar) - 1) / (10 * epsilon * stratPar)
    else:
        tmpLamb = 2 * epsilon
        mistake = 1
        while mistake > 0.01:
            lamb = 2 * epsilon * (1 - 16 * epsilon * tmpLamb * stratPar) ** 0.25
            mistake = abs(tmpLamb - lamb) / tmpLamb
            tmpLamb = lamb
    return lamb


def profileFunc_momentum_Calc(lambInv, epsilon, stratPar):
    '''calculate profile dimensionless function for momentum
        epsilon - main fitting parameter of the model
        stratPar - stratification parameter, mu
        lambInv = 1 / lamb
    '''
    if stratPar > 0:
        return - 5 * epsilon * stratPar / lambInv
    else:
        X = (1 - 16 * epsilon * stratPar / lambInv) ** 0.25
        return 2 * math.log((1 + X) / 2) + \
               math.log((1 + X ** 2) / 2) - 2 * math.atan(X) + math.pi / 2


def profileFunc_heat_Calc(lambInv, epsilon, stratPar):
    '''calculate profile dimensionless function for temperature
        epsilon - main fitting parameter of the model
        stratPar - stratification parameter, mu
        lambInv = 1 / lamb
    '''
    if stratPar > 0:
        return - 5 * epsilon * stratPar / lambInv
    else:
        X = (1 - 16 * epsilon * stratPar / lambInv) ** 0.25
        return 2 * math.log((1 + X ** 2) / 2)


def backgroundCalc(epsilon, M, corPar, MS, MODB, windCpx,
                   tempAtm, tempSurf, roughLength_m,
                   roughLength_h, refLevel):
    '''
    calculate parameters of background atmosphere: friction velocity, temperature scale,
    Ekman depth, geostrophic wind and temperature in free atmosphere
    input:
        epsilon - main fitting parameter of the model,
        M - constant,
        corPar - Coriolis parameter,
        MS - "0" if background is land and "1" if sea,
        MODB - ,
        windCpx - wind at reference level (complex),
        tempAtm - temperature of background atmosphere at reference level,
        tempSurf - background surface temperature,
        roughLength_m - roughness length for momentum for background surface,
        roughLength_h - roughness length for heat for background surface,
        refLevel - reference level, m

    output:
        stratPar - stratification parameter, mu,
        EkmanDepth ,
        frVelCpx - friction velocity for background atmosphere< complex number,
        tempScale - temperature scale for background atmosphere,
        roughLength_m - roughness length for momentum for background surface,
        roughLength_heat - roughness length for heat for background surface,
        tempTop - temperature in free atmosphere,
        gam - ,
        tempSurf - background surface temperature,
        geostrWindCpx - geostrophic wind, complex number
    '''

    CHARNOCKCONST = 0.018
    GRAVITYACC = 9.8
    VONKARMANCONST = 0.41
    BETA = GRAVITYACC / 300.
    VA = 14e-6

    #Iteration for the drag coefficient

    it = 0
    mistake = 1
    DragCoeff = 0.5e-3
    tmpFrVelCpx = math.sqrt(DragCoeff) * windCpx
    stratPar = VONKARMANCONST ** 2. * BETA * (tempAtm - tempSurf) / (corPar * abs(windCpx))
    XX = []
    YY = []
    while mistake > 0.0001:
        it += 1
        tmpFrVel = abs(tmpFrVelCpx)
        roughLength_m = roughLength_m * (1 - MS) + \
                        (CHARNOCKCONST * tmpFrVel ** 2. / GRAVITYACC +
                         0.14 * VA / tmpFrVel) * MS
        lambInv = 1. / lambdaCalc(epsilon, stratPar)
        EkmanDepth = VONKARMANCONST * tmpFrVel / (corPar * lambInv)
        SBLHeight = epsilon * EkmanDepth
        dH = refLevel / EkmanDepth
        dH1 = min(dH, M)
        if dH <= epsilon:
            stratPar1 = dH / epsilon * stratPar
            profileFunc_m = profileFunc_momentum_Calc(lambInv, epsilon, stratPar1)
            profileFunc_h = profileFunc_heat_Calc(lambInv, epsilon, stratPar1)
            frVelCpx = VONKARMANCONST * windCpx /\
                       (math.log(tmpFrVel / (corPar * roughLength_m)) +
                        math.log(VONKARMANCONST * dH / lambInv) - profileFunc_m)
            tempScale = VONKARMANCONST * (tempAtm - tempSurf) /\
                        (math.log(tmpFrVel / (corPar * roughLength_h)) +
                         math.log(VONKARMANCONST * dH / lambInv) - profileFunc_h)
        else:
            profileFunc_m = profileFunc_momentum_Calc(lambInv, epsilon, stratPar)
            profileFunc_h = profileFunc_heat_Calc(lambInv, epsilon, stratPar)
            Coeff = cmath.tanh((1 + 1j) * (M - epsilon))
            Coeff *= (1 - cmath.sinh((1 + 1j) * (M - dH1)) / cmath.sinh((1 + 1j) * (M - epsilon)))
            A = lambInv * Coeff
            B = - A + profileFunc_m  - math.log(VONKARMANCONST * epsilon / lambInv)
            C = -2. * lambInv * (dH - epsilon) + profileFunc_h -\
                math.log(VONKARMANCONST * epsilon / lambInv)
            frVelCpx = VONKARMANCONST * windCpx / \
                       (math.log(tmpFrVel / (corPar * roughLength_m))- B - 1j * A)
            tempScale = VONKARMANCONST * (tempAtm - tempSurf) /\
                                      (math.log(tmpFrVel / (corPar * roughLength_h)) - C)
        mistake = abs(tmpFrVelCpx - frVelCpx) / abs(tmpFrVelCpx)
        tmpFrVelCpx = 0.3 * frVelCpx + 0.7 * tmpFrVelCpx
        tmpFrVel = abs(tmpFrVelCpx)
        stratPar = VONKARMANCONST ** 2 * BETA * tempScale / (corPar * tmpFrVel)
        if it > 100:
            print 'error in backgroundCalc'
            break
    A = lambInv * cmath.tanh((1 + 1j) * (M - epsilon))
    B = - A + profileFunc_momentum_Calc(lambInv, epsilon, stratPar) - math.log(VONKARMANCONST * epsilon / lambInv)
    C = -2 * lambInv * (M - epsilon) + profileFunc_heat_Calc(lambInv, epsilon, stratPar) -\
        math.log(VONKARMANCONST * epsilon / lambInv)
    geostrWindCpx = tmpFrVelCpx / VONKARMANCONST * \
                             (math.log(tmpFrVel / (corPar * roughLength_m))
                              - B - 1j * A)
    tempTop = tempSurf + tempScale / VONKARMANCONST * (math.log(tmpFrVel / (corPar * roughLength_h)) - C)
    gam = 2 * tempScale * abs(frVelCpx) / EkmanDepth ** 2 / corPar
    RES = [stratPar, EkmanDepth, frVelCpx, tempScale, roughLength_m,
           roughLength_h, tempTop, gam, tempSurf, geostrWindCpx]
    return RES


def IBLCalc(epsilon, M, corPar, MS, MODB, stratPar_back, EkmanDepth_back, frVelCpx_back,
            tempScale_back, roughLength_m_back, roughLength_h_back, tempTop, gam, tempSurf_back,
            geostrWindCpx, d, tempSurf, roughLength_h, roughLength_m, stratPar, frVelCpx, tempScale):
    CHARNOCKCONST = 0.018
    GRAVITYACC = 9.8
    VONKARMANCONST = 0.41
    BETA = GRAVITYACC / 300
    VA = 14e-6

    SBLHeight_back = epsilon * EkmanDepth_back
    K_back = EkmanDepth_back ** 2 * corPar / 2
    K = K_back
    geostrWind = abs(geostrWindCpx)
    angle = cmath.phase(geostrWindCpx)
    d -= epsilon
    dm = M - epsilon
    tmpUWind = geostrWind * math.cos(angle)

    # iteration for the drag coefficient
    it = 0
    mistake = 1
    tmpFrVelCpx = frVelCpx
    tmpTempScale = tempScale
    gama = max(0, gam)
    alpha = 1
    epst = 0
    lambInv_back = 1 / lambdaCalc(epsilon, stratPar_back)
    while mistake > 0.0001:
        it += 1
        tmpFrVelCpx = 0.3 * frVelCpx + 0.7 * tmpFrVelCpx
        tmpFrVel = abs(tmpFrVelCpx)
        tmpTempScale = 0.3 * tempScale + 0.7 * tmpTempScale
        heatFlux = - tmpFrVel * tmpTempScale
        stratPar = VONKARMANCONST ** 2 * BETA * tmpTempScale / (corPar * tmpFrVel)
        roughLength_m = roughLength_m * MS + \
                        (CHARNOCKCONST * tmpFrVel ** 2 / GRAVITYACC + 0.14 * VA / tmpFrVel) * (1 - MS)
        lambInv = 1 / lambdaCalc(epsilon, stratPar)
        profileFunc_m = profileFunc_momentum_Calc(lambInv, epsilon, stratPar)
        profileFunc_h = profileFunc_heat_Calc(lambInv, epsilon, stratPar)
        EkmanDepth = VONKARMANCONST * tmpFrVel / (corPar * lambInv)
        SBLHeight = epsilon * EkmanDepth
        K = EkmanDepth ** 2 * corPar / 2
        delta = EkmanDepth * (d + epsilon)
        if delta <= SBLHeight_back:
            tmpStratPar_back = stratPar_back * delta / SBLHeight_back
            profileFunc_m_back = profileFunc_momentum_Calc(lambInv_back, epsilon, tmpStratPar_back)
            profileFunc_h_back = profileFunc_heat_Calc(lambInv_back, epsilon, tmpStratPar_back)
            windCpx_back = frVelCpx_back / VONKARMANCONST * (math.log(delta / roughLength_m_back) - profileFunc_m_back)
            temp = tempSurf_back + tempScale_back / VONKARMANCONST *\
                   (math.log(delta / roughLength_h_back) - profileFunc_h_back)
        else:
            ksi = (delta - SBLHeight_back) / (EkmanDepth_back * dm)
            tmpKsi = ksi
            ksi = min(1, ksi)
            Coeff = (1 -1j) * cmath.sinh((1 + 1j) * dm * (1 - ksi)) / cmath.cosh((1 + 1j) * dm)
            windCpx_back = geostrWindCpx - lambInv_back * frVelCpx_back / VONKARMANCONST * Coeff
            temp = tempTop - 2 * lambInv_back * tempScale_back / VONKARMANCONST * dm * (1 - tmpKsi)
        if gam <= 0:
            alg = 1
            epst = 0
        else:
            alg = (abs(heatFlux) + abs(gama * K_back)) / (abs(heatFlux) + abs(gama * K))
            alg = min(1, alg)
            epst = max(0, 1 / 4 * heatFlux / (K * gama * alg))
        temp -= epst * gama * delta
        alpha = alg * (1 - (delta / (M * EkmanDepth)) ** 4)
        alpha = max(0, alpha)
        Coeff = cmath.tanh((1 + 1j) * d)
        A = lambInv * Coeff
        B = - A + profileFunc_m  - math.log(VONKARMANCONST * epsilon / lambInv)
        C = -2 * lambInv * d + profileFunc_h - math.log(VONKARMANCONST * epsilon / lambInv)

        F = 1j * VONKARMANCONST / tmpFrVelCpx * BETA / (corPar * tmpUWind) * d ** 2 / (d ** 2 - 1j * alpha)
        F *= heatFlux
        F *= MODB
        frVelCpx = VONKARMANCONST * geostrWindCpx / (math.log(tmpFrVel / (corPar * roughLength_m)) - B - 1j * A + F)
        frVelCpx *= (1 + (windCpx_back / geostrWindCpx -1) / cmath.cosh((1 + 1j) * d))
        tempScale = VONKARMANCONST * (temp - tempSurf) / (math.log(tmpFrVel / (corPar * roughLength_h)) - C)
        mistake = abs(tmpFrVelCpx - frVelCpx) / abs(tmpFrVelCpx)
        if it > 100:
            print 'error in IBLCalc'
            break
    RES = [stratPar, EkmanDepth, frVelCpx, tempScale, temp, roughLength_m, windCpx_back, alpha, epst]
    return RES


def IBLLogCalc(epsilon, M, corPar, MS, MODB, stratPar_back, EkmanDepth_back, frVelCpx_back,
            tempScale_back, roughLength_m_back, roughLength_h_back, tempTop, gam, tempSurf_back,
            geostrWindCpx, d, tempSurf, roughLength_h, roughLength_m, stratPar, frVelCpx, tempScale, MOD=None):

    CHARNOCKCONST = 0.018
    GRAVITYACC = 9.8
    VONKARMANCONST = 0.41
    BETA = GRAVITYACC / 300
    VA = 14e-6

    SBLHeight_back = epsilon * EkmanDepth_back
    K_back = EkmanDepth_back ** 2 * corPar / 2
    K = K_back
    geostrWind = abs(geostrWindCpx)
    angle = cmath.phase(geostrWindCpx)

    # iteration for the drag coefficient
    it = 0
    mistake = 1
    tmpFrVelCpx = frVelCpx
    tmpTempScale = tempScale
    K_back = EkmanDepth_back ** 2 * corPar / 2
    gama = max(0, gam)
    lambInv_back = 1 / lambdaCalc(epsilon, stratPar_back)
    while mistake > 0.0001:
        it += 1
        tmpFrVelCpx = 0.3 * frVelCpx + 0.7 * tmpFrVelCpx
        tmpFrVel = abs(tmpFrVelCpx)
        tmpTempScale = 0.3 * tempScale + 0.7 * tmpTempScale
        heatFlux = - tmpFrVel * tempScale
        stratPar = VONKARMANCONST ** 2 * BETA * tmpTempScale / (corPar * tmpFrVel)
        roughLength_m = roughLength_m * MS + \
                        (CHARNOCKCONST * tmpFrVel ** 2 / GRAVITYACC + 0.14 * VA / tmpFrVel) * (1 - MS)
        lambInv = 1 / lambdaCalc(epsilon, stratPar)
        if MOD != 'initial':
            tmpStratPar = stratPar * d / epsilon
        else:
            tmpStratPar = stratPar * d / EkmanDepth_back / epsilon
        profileFunc_m = profileFunc_momentum_Calc(lambInv, epsilon, tmpStratPar)
        profileFunc_h = profileFunc_heat_Calc(lambInv, epsilon, tmpStratPar)
        EkmanDepth = VONKARMANCONST * tmpFrVel / (corPar * lambInv)
        SBLHeight = epsilon * EkmanDepth
        K = EkmanDepth ** 2 * corPar / 2
        if MOD != 'initial':
            delta = EkmanDepth * d
        else:
            delta = d
        if delta <= SBLHeight_back:
            tmpStratPar_back = stratPar_back * delta / SBLHeight_back
            profileFunc_m_back = profileFunc_momentum_Calc(lambInv_back, epsilon, tmpStratPar_back)
            profileFunc_h_back = profileFunc_heat_Calc(lambInv_back, epsilon, tmpStratPar_back)
            windCpx_back = frVelCpx_back / VONKARMANCONST * (math.log(delta / roughLength_m_back) - profileFunc_m_back)
            temp = tempSurf_back + tempScale_back / VONKARMANCONST *\
                   (math.log(delta / roughLength_h_back) - profileFunc_h_back)
        else:
            dm = M - epsilon
            ksi = (delta - SBLHeight_back) / (EkmanDepth_back * dm)
            tmpKsi = ksi
            ksi = min(1, ksi)
            Coeff = (1 -1j) * cmath.sinh((1 + 1j) * dm * (1 - ksi)) / cmath.cosh((1 + 1j) * dm)
            windCpx_back = geostrWindCpx - lambInv_back * frVelCpx_back / VONKARMANCONST * Coeff
            temp = tempTop - 2 * lambInv_back * tempScale_back / VONKARMANCONST * dm * (1 - tmpKsi)
        if gam <= 0:
            alg = 1
            epst = 0
        else:
            alg = (abs(heatFlux) + abs(gama * K_back)) / (abs(heatFlux) + abs(gama * K))
            alg = min(1, alg)
            epst = max(0, 1 / 4 * heatFlux / (K * gama * alg))
        temp -= epst * gama * delta
        frVelCpx = VONKARMANCONST * windCpx_back / (math.log(delta / roughLength_m) - profileFunc_m)
        tempScale = VONKARMANCONST * (temp - tempSurf) / (math.log(delta / roughLength_h) - profileFunc_h)
        mistake = abs(tmpFrVelCpx - frVelCpx) / abs(tmpFrVelCpx)
        if it > 100:
            print 'error in IBLCalc'
            alpha = alg * (1 - (delta / (M * EkmanDepth)) ** 4)
            break
    alpha = alg * (1 - (delta / (M * EkmanDepth)) ** 4)
    RES = [stratPar, EkmanDepth, frVelCpx, tempScale, temp, roughLength_m, windCpx_back, alpha, epst]
    return RES

def IBLModelProcessing():
    CHARNOCKCONST = 0.018
    GRAVITYACC = 9.8
    OMEGA = 7.3e-5
    VONKARMANCONST = 0.41
    BETA = GRAVITYACC / 300
    VA = 14e-6
    lat = 56.


    corPar = 2. * math.sin(math.radians(lat)) * OMEGA
    epsilon = 0.15
    M = 2.05
    MODB = 0
    ##
    MODP = 1
    refLevel = 10
    windSpeed = 5
    angle = 30
    MS = 0
    tempSurf_back = 20
    tempAtm_back = 15
    tempSurf = 18
    roughLength_m_back = 0.1
    roughLength_m = 0
    ##
    angle = math.radians(angle)
    windCpx = windSpeed * cmath.exp(1j * angle)
    if MS == 0:
        roughLength_h_back = roughLength_m_back * math.exp(-2. / 3.)
        roughLength_m = 0
        roughLength_h = 1e-5
    else:
        roughLength_m_back = 0
        roughLength_h_back = 1e-5
        roughLength_h = roughLength_m * math.exp(-2. / 3.)
    dTemp = tempSurf - tempSurf_back
    # calculating parameters of background atmosphere
    inputArray = (epsilon, M, corPar, MS, MODB, windCpx,
                   tempAtm_back, tempSurf_back, roughLength_m_back,
                   roughLength_h_back, refLevel)
    output = backgroundCalc(*inputArray)
    (stratPar_back, EkmanDepth_back, frVelCpx_back, tempScale_back, roughLength_m_back,
    roughLength_h_back, tempTop, gam, tempSurf_back, geostrWindCpx) = output
    geostrWind = abs(geostrWindCpx)
    angle = cmath.phase(geostrWindCpx)

    LS = geostrWind / corPar

    # initiate dparmin
    stratPar = VONKARMANCONST ** 2. * BETA * (tempAtm_back - tempSurf) / (corPar * windSpeed)
    frVelCpx = math.sqrt(1e-3) * windCpx
    tempScale = math.sqrt(1e-3) * (tempAtm_back - tempSurf)
    delta = 10. * roughLength_m_back
    inputArray = (epsilon, M, corPar, MS, MODB, stratPar_back, EkmanDepth_back, frVelCpx_back,
            tempScale_back, roughLength_m_back, roughLength_h_back, tempTop, gam, tempSurf_back,
            geostrWindCpx, delta, tempSurf, roughLength_h, roughLength_m, stratPar, frVelCpx, tempScale, 'initial')
    output = IBLLogCalc(*inputArray)
    (stratPar, EkmanDepth, frVelCpx, tempScale, temp, roughLength_m, windCpx_back, alpha, epst) = output
    dparmin = delta / EkmanDepth
    ND = 15
    dparmax = M - 0.05
    dlnd = (math.log(dparmax) - math.log(dparmin)) / (ND - 1)
    dpar = [dparmin * math.exp(x * dlnd) for x in xrange(0, ND)]
    # loop over d
    F = []
    coeff_G = []
    coeff_TR = []
    dAngle = []
    EkmanDepth_array = []
    frVelCpx_array = []
    tempScale_array = []
    delta_array = []
    SBLHeight = []
    roughLength_m_array = []
    stratPar_array = []
    windCpx_back_array = []
    temp_array = []
    l_array = []
    alpha_array = []
    epst_array = []
    for d in dpar:
        if d <= epsilon:
            inputArray = (epsilon, M, corPar, MS, MODB, stratPar_back, EkmanDepth_back, frVelCpx_back,
                    tempScale_back, roughLength_m_back, roughLength_h_back, tempTop, gam, tempSurf_back,
                    geostrWindCpx, d, tempSurf, roughLength_h, roughLength_m, stratPar, frVelCpx, tempScale)
            output = IBLLogCalc(*inputArray)
            (stratPar, EkmanDepth, frVelCpx, tempScale, temp, roughLength_m, windCpx_back, alpha, epst) = output
            c1 = epsilon / d * lambdaCalc(epsilon, stratPar) / lambdaCalc(epsilon, stratPar * d / epsilon)
            c2 = abs(windCpx_back) * math.cos(cmath.phase(windCpx_back)) / geostrWind
            F.append(c1 * c2 / EkmanDepth * d / alpha)
        else:
            inputArray = (epsilon, M, corPar, MS, MODB, stratPar_back, EkmanDepth_back, frVelCpx_back,
                    tempScale_back, roughLength_m_back, roughLength_h_back, tempTop, gam, tempSurf_back,
                    geostrWindCpx, d, tempSurf, roughLength_h, roughLength_m, stratPar, frVelCpx, tempScale)
            output = IBLCalc(*inputArray)
            (stratPar, EkmanDepth, frVelCpx, tempScale, temp, roughLength_m, windCpx_back, alpha, epst) = output
            lambInv = 1 / lambdaCalc(epsilon, stratPar)
            profileFunc_m = profileFunc_momentum_Calc(lambInv, epsilon, stratPar)
            c1 = frVelCpx / VONKARMANCONST * (math.log(epsilon * EkmanDepth / roughLength_m) - profileFunc_m)
            c2 = (c1 + windCpx_back) / 2.
            c3 = abs(c2) * math.cos(cmath.phase(c2)) / geostrWind
            F.append(c3 / EkmanDepth * d / alpha)
        coeff_G.append((abs(frVelCpx) / abs(geostrWindCpx)) ** 2)
        if (temp - tempSurf) == 0:
            coeff_TR.append(None)
        else:
            coeff_TR.append((abs(frVelCpx) / abs(geostrWindCpx) * tempScale / (temp - tempSurf)))
        dAngle.append((cmath.phase(frVelCpx) - cmath.phase(geostrWindCpx)) * 180. / math.pi)
        EkmanDepth_array.append(EkmanDepth)
        frVelCpx_array.append(frVelCpx)
        tempScale_array.append(tempScale)
        delta_array.append(d * EkmanDepth)
        SBLHeight.append(epsilon * EkmanDepth)
        roughLength_m_array.append(roughLength_m)
        stratPar_array.append(stratPar)
        windCpx_back_array.append(windCpx_back)
        temp_array.append(temp)
        if stratPar == 0:
            l_array.append(None)
        else:
            l_array.append(VONKARMANCONST * abs(frVelCpx) / (corPar * stratPar))
        alpha_array.append(alpha)
        epst_array.append(epst)
    # calculation of FETCH X for each delta, x is dimentionless fetch normalized on LS= geostrWind / corPar
    X = [delta_array[0] * math.log(delta_array[0] / roughLength_m_array[0]) / (2. * VONKARMANCONST ** 2 * LS)]
    for i in xrange(1, len(delta_array)):
        X.append((F[i] + F[i - 1]) / 2. * abs(delta_array[i] - delta_array[i - 1]) + X[i - 1])
    X = [x * LS for x in X]
    if MODP == 1:
        inputArray = [EkmanDepth_array, EkmanDepth_back, delta_array, frVelCpx_array, frVelCpx_back,
                 stratPar_array, stratPar_back, SBLHeight, roughLength_m_array, roughLength_m_back,
                 roughLength_h, roughLength_h_back, tempSurf, tempTop, temp_array, tempSurf_back,
                 tempScale_array, tempScale_back, windCpx_back_array, geostrWindCpx, dpar,
                 alpha_array, corPar, angle, epsilon, M, MODB]
        wind_prof, temp_prof, z_array, wind_prof_back, temp_prof_back = profileCalc(*inputArray)
        plotProfiles(wind_prof, wind_prof_back, z_array, 'Wind profiles')
        plotProfiles(temp_prof, temp_prof_back, z_array, 'Temperature profiles')
    #    plt.show()
    plotFetches(X, delta_array, SBLHeight, 'SBL and IBL height as function of fetch')
    plot4(X, coeff_G, coeff_TR, dAngle, delta_array, SBLHeight)
    plt.show()


def profileCalc(EkmanDepth_array, EkmanDepth_back, delta_array, frVelCpx_array, frVelCpx_back,
                stratPar_array, stratPar_back, SBLHeight_array, roughLength_m_array, roughLength_m_back,
                roughLength_h, roughLength_h_back, tempSurf, tempTop, temp_array, tempSurf_back,
                tempScale_array, tempScale_back, windCpx_back_array, geostrWindCpx, dpar,
                alpha_array, corPar, angle, epsilon, M, MODB):

    CHARNOCKCONST = 0.018
    GRAVITYACC = 9.8
    OMEGA = 7.29e-5
    VONKARMANCONST = 0.41
    BETA = GRAVITYACC / 300
    N = 50

    zmax = max(max([x * 2 for x in EkmanDepth_array]), EkmanDepth_back)
    hmax = zmax
    hmin = 1
    dlnz = math.log(hmax / hmin) / (N - 1)
    J = 0
    wind_prof = [[] for x in xrange(0, len(delta_array), 2)]
    temp_prof = [[] for x in xrange(0, len(delta_array), 2)]
    z_array = [[] for x in xrange(0, len(delta_array), 2)]
    for i in xrange(0, len(delta_array), 2):
        d = dpar[i] - epsilon
        frVel = abs(frVelCpx_array[i])
        lambInv = VONKARMANCONST * frVel / (EkmanDepth_array[i] * corPar)
        profileFunc_m = profileFunc_momentum_Calc(lambInv, epsilon, stratPar_array[i])
        u_h = frVelCpx_array[i] / VONKARMANCONST * (math.log(SBLHeight_array[i] / roughLength_m_array[i]) - profileFunc_m)
        wind_prof_back = []
        temp_prof_back = []
        z = []
        for indx in xrange(N):
            z.append(hmin * math.exp(dlnz * indx))
            if z[indx] <= epsilon * EkmanDepth_back:
                stratPar = stratPar_back * z[indx] / (epsilon * EkmanDepth_back)
                tmpLambInv = VONKARMANCONST * abs(frVelCpx_back) / (EkmanDepth_back * corPar)
                profileFunc_m = profileFunc_momentum_Calc(tmpLambInv, epsilon, stratPar)
                profileFunc_h = profileFunc_heat_Calc(tmpLambInv, epsilon, stratPar)
                wind_prof_back.append(frVelCpx_back / VONKARMANCONST *
                                 (math.log(z[indx] / roughLength_m_back) - profileFunc_m))
                temp_prof_back.append(tempSurf_back + tempScale_back / VONKARMANCONST *
                                 (math.log(z[indx] / roughLength_h_back) - profileFunc_h))
            else:
                dmax = M - epsilon
                ksi = (z[indx] - epsilon * EkmanDepth_back) / (EkmanDepth_back * dmax)
                tmpKsi = ksi
                ksi = min(1, ksi)
                coeff = (1 - 1j) * cmath.sinh((1 + 1j) * dmax * (1 - ksi)) / cmath.cosh((1 + 1j) * dmax)
                wind_prof_back.append(geostrWindCpx - 1. / lambdaCalc(epsilon, stratPar_back) *
                                                frVelCpx_back / VONKARMANCONST * coeff)
                temp_prof_back.append(tempTop - 2. / lambdaCalc(epsilon, stratPar_back) *
                                           tempScale_back / VONKARMANCONST * dmax * (1 - tmpKsi))
            if z[indx] <= delta_array[i]:
                if z[indx] <= SBLHeight_array[i]:
                    stratPar = stratPar_array[i] * z[indx] / SBLHeight_array[i]
                    profileFunc_m = profileFunc_momentum_Calc(lambInv, epsilon, stratPar)
                    profileFunc_h = profileFunc_heat_Calc(lambInv, epsilon, stratPar)
                    wind_prof[J].append(frVelCpx_array[i] / VONKARMANCONST *
                                     (math.log(z[indx] / roughLength_m_array[i]) - profileFunc_m))
                    temp_prof[J].append(tempSurf + tempScale_array[i] / VONKARMANCONST *
                                     (math.log(z[indx] / roughLength_h) - profileFunc_h))
                else:
                    ksi = (z[indx] - SBLHeight_array[i]) / (delta_array[i] - SBLHeight_array[i])
                    tmpKsi = ksi
                    tmpKsi = ksi
                    ksi = min(1, ksi)
                    heatFlux = - frVel * tempScale_array[i]
                    alpha = alpha_array[i]
                    coeff1 = (1 - 1j) * cmath.sinh((1 + 1j) * d * (1 - ksi)) / cmath.cosh((1 + 1j) * d)
                    coeff2 = heatFlux * (1 - ksi)
                    tmpCoeff2 = 0
                    coeff3 = heatFlux * (1 - ksi)
                    coeff4 = BETA * d ** 2 / (corPar * abs(geostrWindCpx) * math.cos(angle) * (alpha + 1j * d ** 2))
                    coeff4 *= (coeff2 - tmpCoeff2)
                    coeff4 *= MODB
                    coeff5 = (windCpx_back_array[i] - geostrWindCpx) * \
                             cmath.cosh((1 + 1j) * d * ksi) / cmath.cosh((1 + 1j) * d)
                    wind_prof[J].append(geostrWindCpx - lambInv * frVelCpx_array[i] / VONKARMANCONST * coeff1\
                                   + coeff4 +coeff5)
                    temp_prof[J].append(temp_array[i] + 2. * lambInv * d * coeff3 / (frVel * VONKARMANCONST))
            else:
                wind_prof[J].append(wind_prof_back[indx])
                temp_prof[J].append(temp_prof_back[indx])
            z_array[J].append(z[indx])
        J += 1
    return wind_prof, temp_prof, z_array, wind_prof_back, temp_prof_back


def plotProfiles(x, x2, y, title=None):
    fig = plt.figure(facecolor='w')
    ax = fig.add_subplot(1, 1, 1)
    for i in xrange(len(x)):
        plt.plot([abs(xx) for xx in x[i]], y[i])
    plt.plot([abs(xx) for xx in x2], y[i], '--w', linewidth=2.0)
    ax.set_yscale('log')
    ax.grid(True, which='both')
    if title:
        plt.title(title)
    return


def plotFetches(x, y1, y2, title=None):
    fig = plt.figure(facecolor='w')
    ax = fig.add_subplot(1, 1, 1)
    plt.plot(x, y1)
    plt.plot(x, y2)
    #plt.plot([abs(xx) for xx in x2], y[i], '--w', linewidth=2.0)
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.grid(True, which='both')
    if title:
        plt.title(title)
    return


def plot4(x, y1, y2, y3, y4, y5):
    fig = plt.figure(facecolor='w')
    ax = fig.add_subplot(2, 2, 1)
    plt.plot(x, y1)
    plt.title('geostrophic drag coefficient vs. fetch')
    ax.set_xscale('log')
    ax.grid(True, which='both')
    ax = fig.add_subplot(2, 2, 2)
    plt.plot(x, y2)
    plt.title('termal resistance coefficient vs. fetch')
    ax.set_xscale('log')
    ax.grid(True, which='both')
    ax = fig.add_subplot(2, 2, 3)
    plt.plot(x, y3)
    plt.title('wind direction turn vs. fetch')
    ax.set_xscale('log')
    ax.grid(True, which='both')
    ax = fig.add_subplot(2, 2, 4)
    plt.plot(x, y4)
    plt.plot(x, y5)
    plt.title('SBL and IBL height as function of fetch')
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.grid(True, which='both')

    return



def main():
    IBLModelProcessing()

if __name__ == "__main__":
    main()