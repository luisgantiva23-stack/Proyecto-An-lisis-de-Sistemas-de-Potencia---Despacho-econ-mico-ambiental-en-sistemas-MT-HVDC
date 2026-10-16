using JuMP
using Ipopt
using Printf

# ==========================================
# 1. DEFINICIÓN DE PARÁMETROS BASE
# ==========================================
S_base = 1000.0 # Potencia base en MVA (o MW)
V_base = 400.0  # Tensión base en kV
Z_base = (V_base^2) / S_base # Impedancia base en Ohms (160 Ω)
I_base = S_base / V_base     # Corriente base en kA (2.5 kA)

num_nodos = 6
num_periodos = 4
num_gen = 3
nodo_slack = 2

# ==========================================
# 2. DATOS DEL SISTEMA
# ==========================================
# Datos de Generadores: [c ($), b ($/MWh), a ($/MW^2), Pmin (MW), Pmax (MW)]
gen_costos = [
    100.0  20.0  0.10   50.0 1500.0;
    100.0  15.0  0.12  100.0 2000.0;
    200.0  18.0  0.04  140.0 1800.0
]

# Datos de Emisiones: [gamma (kg), beta (kg/MWh), alpha (kg/MW^2)]
gen_emisiones = [
    4.091  -5.543  0.06490;
    2.543  -6.047  0.05638;
    4.258  -5.094  0.04586
]

# Topología de Líneas: [Desde, Hacia, Resistencia (Ohms)]
lineas_data = [
    1 5 5.70;
    5 3 2.28;
    5 4 1.71;
    1 3 2.28;
    3 6 4.75;
    1 2 1.90;
    2 6 1.90
]
num_lineas = size(lineas_data, 1)

# Límite térmico de líneas
I_max_kA = 4.75 
I_max_pu = I_max_kA / I_base

# Demandas por nodo (4, 5 y 6) en los 4 periodos (en MW)
Pd_MW = zeros(num_nodos, num_periodos)
Pd_MW[4, :] = [1125.0, 1200.0, 1500.0, 750.0]
Pd_MW[5, :] = [1250.0, 1000.0, 750.0,  875.0]
Pd_MW[6, :] = [665.0,  855.0,  760.0,  950.0]

# Conversión a p.u.
Pd_pu = Pd_MW ./ S_base

# Construcción de la matriz de Conductancias (G) en p.u.
G = zeros(num_nodos, num_nodos)
for i in 1:num_lineas
    from = Int(lineas_data[i, 1])
    to = Int(lineas_data[i, 2])
    R_pu = lineas_data[i, 3] / Z_base
    g_line = 1.0 / R_pu
    
    G[from, to] -= g_line
    G[to, from] -= g_line
    G[from, from] += g_line
    G[to, to] += g_line
end

# ==========================================
# 3. FUNCIÓN GENERAL DEL DESPACHO PONDERADO
# ==========================================
function resolver_despacho(w::Float64, nombre_escenario::String)
    modelo = Model(Ipopt.Optimizer)
    set_silent(modelo) # Ocultar logs iterativos del solver
    
    # Variables de decisión
    @variable(modelo, 0.90 <= v[1:num_nodos, 1:num_periodos] <= 1.10)
    @variable(modelo, p_g[1:num_gen, 1:num_periodos] >= 0)
    
    # Fijar voltaje del Nodo Slack (con force=true)
    for h in 1:num_periodos
        fix(v[nodo_slack, h], 1.0; force = true)
    end
    
    # Límites de capacidad técnica de los generadores
    for i in 1:num_gen
        for h in 1:num_periodos
            set_lower_bound(p_g[i, h], gen_costos[i, 4] / S_base)
            set_upper_bound(p_g[i, h], gen_costos[i, 5] / S_base)
        end
    end
    
    # FUNCIÓN OBJETIVO UNIFICADA
    @objective(modelo, Min, 
        w * sum(
            gen_costos[i, 3] * (p_g[i, h] * S_base)^2 +
            gen_costos[i, 2] * (p_g[i, h] * S_base) +
            gen_costos[i, 1]
            for i in 1:num_gen, h in 1:num_periodos
        )
        +
        (1.0 - w) * sum(
            gen_emisiones[i, 3] * (p_g[i, h] * S_base)^2 +
            gen_emisiones[i, 2] * (p_g[i, h] * S_base) +
            gen_emisiones[i, 1]
            for i in 1:num_gen, h in 1:num_periodos
        )
    )
    
    # RESTRICCIONES DE RED
    for h in 1:num_periodos
        # Balance de Potencia 
        for k in 1:num_nodos
            P_gen_k = (k <= 3) ? p_g[k, h] : 0.0
            @NLconstraint(modelo, P_gen_k - Pd_pu[k, h] == v[k, h] * sum(G[k, m] * v[m, h] for m in 1:num_nodos))
        end
        
        # Límites térmicos de corriente en las líneas
        for i in 1:num_lineas
            from = Int(lineas_data[i, 1])
            to = Int(lineas_data[i, 2])
            R_pu = lineas_data[i, 3] / Z_base
            @constraint(modelo, (v[from, h] - v[to, h]) / R_pu <= I_max_pu)
            @constraint(modelo, (v[from, h] - v[to, h]) / R_pu >= -I_max_pu)
        end
    end
    
    # ==========================================
    # 4. EXTRACCIÓN E IMPRESIÓN DE RESULTADOS
    # ==========================================
    optimize!(modelo)
    
    if termination_status(modelo) == LOCALLY_SOLVED || termination_status(modelo) == OPTIMAL
        P_opt_MW = value.(p_g) .* S_base
        V_opt_pu = value.(v)
        
        costo_total = sum(
            gen_costos[i, 3] * (P_opt_MW[i, h])^2 +
            gen_costos[i, 2] * (P_opt_MW[i, h]) +
            gen_costos[i, 1]
            for i in 1:num_gen, h in 1:num_periodos
        )
        
        emisiones_totales = sum(
            gen_emisiones[i, 3] * (P_opt_MW[i, h])^2 +
            gen_emisiones[i, 2] * (P_opt_MW[i, h]) +
            gen_emisiones[i, 1]
            for i in 1:num_gen, h in 1:num_periodos
        )

        # ------------------------------------
        # BLOQUE DE IMPRESIÓN
        # ------------------------------------
        println("="^65)
        println("REPORTE: $nombre_escenario (w = $w)")
        println("="^65)
        
        @printf("Costo Operativo Total   : \$ %.2f USD\n", costo_total)
        @printf("Emisiones Equivalentes  : %.2f kg CO2\n\n", emisiones_totales)
        
        println("1. Despacho de Potencia por Generador (MW):")
        println("Generador\th1\t\th2\t\th3\t\th4")
        for i in 1:num_gen
            @printf("P%d\t\t%.4f\t%.4f\t%.4f\t%.4f\n", i, P_opt_MW[i, 1], P_opt_MW[i, 2], P_opt_MW[i, 3], P_opt_MW[i, 4])
        end
        
        println("\n2. Perfil de Tensiones por Nodo (%):")
        println("Nodo\th1\t\th2\t\th3\t\th4")
        for k in 1:num_nodos
            @printf("%d\t%.4f\t\t%.4f\t\t%.4f\t\t%.4f\n", k, V_opt_pu[k, 1]*100, V_opt_pu[k, 2]*100, V_opt_pu[k, 3]*100, V_opt_pu[k, 4]*100)
        end
        
        println("\n3. Cargabilidad de Líneas (%% respecto a $(I_max_kA) kA):")
        println("Línea\th1\t\th2\t\th3\t\th4")
        for i in 1:num_lineas
            from = Int(lineas_data[i, 1])
            to = Int(lineas_data[i, 2])
            R_pu = lineas_data[i, 3] / Z_base
            
            @printf("%d\t", i)
            for h in 1:num_periodos
                I_flujo_pu = abs((V_opt_pu[from, h] - V_opt_pu[to, h]) / R_pu)
                porcentaje = (I_flujo_pu / I_max_pu) * 100
                @printf("%.2f\t\t", porcentaje)
            end
            println("")
        end
        println("\n")
    else
        println("El modelo no encontró una solución óptima para: $nombre_escenario.")
    end
end

# ==========================================
# 5. EJECUCIÓN DE LOS 3 ESCENARIOS
# ==========================================
resolver_despacho(1.0, "Escenario 1 - Minimización de Costos")
resolver_despacho(0.0, "Escenario 2 - Minimización de Emisiones")
resolver_despacho(0.5, "Escenario 3 - Ponderación Equitativa 50/50")