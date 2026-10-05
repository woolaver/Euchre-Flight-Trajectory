syms a t

func_1 = a*t == 2652;
func_2 = .5*a*t^2 == 18.288;

[a_sol, t_sol] = solve([func_1, func_2], [a, t]);

disp('Required Acceleration:')
disp(double(a_sol))
disp(double(t_sol))

g_force = a_sol/9.8;
disp(double(g_force))