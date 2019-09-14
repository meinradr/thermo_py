from thermo_py import RTD


rtd_1 = RTD()
rtd_2 = RTD(r0=99)

t_values = [-50, -25, 0, 50, 100]
rt_values = [80, 100, 120, 140]

c_rt = RTD.compare_rt(rtd_1, rtd_2, t_values)
c_t = RTD.compare_t(rtd_1, rtd_2, rt_values)

print(c_rt)
print(c_t)




