from iapws import IAPWS97
import numpy as np

# 给定初始参数
Ne = 1000
eta_1 = 0.995   # 0.99-1.0
xfh = 0.9975
xi_d = 0.0105
eta_h_i = 0.8207
eta_l_i = 0.8359
eta_m = 0.985   # 0.98-0.99
eta_ge = 0.985   # 0.98-0.99
theta_h_u = 3
theta_l_u = 2
eta_h = 0.98   # 0.97-0.99
eta_fwp_p = 0.58
eta_fwp_ti = 0.80   # 0.78-0.82
eta_fwp_tm = 0.90
eta_fwp_tg = 0.98
T_sw_1 = 24

# 蒸汽初参数
pc = 15.5   # 15-16
print("1.反应堆冷却剂系统运行压力", format(pc, ".1f"))
Tc_s = IAPWS97(P=pc, x=0).T - 273.15
print("2.冷却剂压力对应的饱和温度", format(Tc_s, ".2f"))
delta_Tsub = 17.5   # 15-20
print("3.反应堆出口冷却剂过冷度", format(delta_Tsub, ".1f"))
Tco = Tc_s - delta_Tsub
print("4.反应堆出口冷却剂温度", format(Tco, ".2f"))
delta_Tc = 35   # 30-40
print("5.反应堆进出口冷却剂温升", format(delta_Tc, ".0f"))
Tci = Tco - delta_Tc
print("6.反应堆进口冷却剂温度", format(Tci, ".2f"))
ps = 6  # 5.0-7.0
print("7.蒸汽发生器饱和蒸汽压力", format(ps, ".1f"))
Ts = IAPWS97(P=ps, x=1).T - 273.15
print("8.蒸汽发生器饱和蒸汽温度", format(Ts, ".2f"))
delta_Tm = (Tco - Tci) / np.log((Tco - Ts) / (Tci - Ts))
print("9.一、二次侧对数平均温差", format(delta_Tm, ".2f"))

# 蒸汽终参数
delta_Tsw = 7   # 6-8
print("10.冷凝器中循环冷却水温升", format(delta_Tsw, ".0f"))
delta_T = 6.5   # 3-10
print("11.冷凝器传热端差", format(delta_T, ".1f"))
Tcd = T_sw_1 + delta_Tsw + delta_T
print("12.冷凝器凝结水饱和温度", format(Tcd, ".2f"))
pcd = IAPWS97(T=Tcd + 273.15, x=0).P
print("13.冷凝器的运行压力", format(pcd, ".4f"))

# 蒸汽中间再热参数
pfh = ps
delta_pfh = 0.05 * pfh     # 0.03-0.07
ph_i = pfh - delta_pfh
print("14.高压缸进口蒸汽压力", format(ph_i, ".2f"))
hfh = IAPWS97(P=pfh, x=xfh).h
hhi = hfh
xh_i = IAPWS97(P=ph_i, h=hfh).x
print("15.高压缸进口蒸汽干度", format(xh_i, ".4f"))
ph_z = 0.13 * ph_i  # 0.12-0.14
print("16.高压缸排汽压力", format(ph_z, ".2f"))
sh_i = IAPWS97(P=ph_i, x=xh_i).s
hih_z = IAPWS97(P=ph_z, s=sh_i).h
hh_z = hfh - eta_h_i * (hfh - hih_z)
xh_z = IAPWS97(P=ph_z, h=hh_z).x
print("17.高压缸排汽干度", format(xh_z, ".4f"))
psp_i = ph_z
print("18.汽水分离器进口蒸汽压力", format(psp_i, ".2f"))
xsp_i = xh_z
print("19.汽水分离器进口蒸汽干度", format(xsp_i, ".4f"))
puw = 0.99 * ph_z  # 汽水分离器出口疏水压力
huw = IAPWS97(P=puw, x=0).h  # 汽水分离器出口疏水比焓

# 第一级再热器
prh1_i = 0.99 * ph_z
print("20.再热蒸汽进口压力", format(prh1_i, ".4f"))
xrh1_i = xsp_i / (1 - 0.98 * (1 - xsp_i))   # 参照大亚湾的蒸汽参数，汽水分离器除去蒸汽中98%的水
print("21.再热蒸汽进口干度", format(xrh1_i, ".4f"))
prh1_hs = ph_z
print("22.加热蒸汽进口压力", format(prh1_hs, ".4f"))
xrh1_hs = xh_z
print("23.加热蒸汽进口干度", format(xrh1_hs, ".4f"))

# 第二级再热器
prh2_i = 0.98 * ph_z    # 再热蒸汽进口压力
prh2_z = 0.97 * ph_z    # 二级再热器出口压力
Trh2_z = Ts - 14    # 二级再热器出口温度 13-15
Trh1_i = IAPWS97(P=prh1_i, x=xrh1_i).T - 273.15
hrh1_i = IAPWS97(P=prh1_i, x=xrh1_i).h
hrh2_z = IAPWS97(P=prh2_z, T=Trh2_z+273.15).h   # 二级再热器出口蒸汽比焓
delta_hrh = (hrh2_z - hrh1_i)/2     # 每级再热器平均焓升
hli = hrh2_z
hrh1_z = (hrh1_i + hrh2_z) / 2  # 一级再热器出口蒸汽比焓
Trh2_i = IAPWS97(P=prh2_i, h=hrh1_z).T - 273.15     # 二级再热器进口蒸汽温度
prh2_hs = ph_i  # 加热蒸汽进口压力
xrh2_hs = xh_i  # 加热蒸汽进口干度

print("24.再热蒸汽进口压力", format(prh2_i, ".4f"))
print("25.再热蒸汽进口温度", format(Trh2_i, ".2f"))
print("26.再热蒸汽出口压力", format(prh2_z, ".4f"))
print("27.再热蒸汽出口温度", format(Trh2_z, ".2f"))
print("28.加热蒸汽进口压力", format(prh2_hs, ".2f"))
print("29.加热蒸汽进口干度", format(xrh2_hs, ".4f"))

# 低压缸
pl_i = prh2_z * 0.98    # 考虑低压缸的进汽损失占再热器出口压力的2%
print("30.进口蒸汽压力", format(pl_i, ".4f"))
Tl_i = IAPWS97(P=pl_i, h=hrh2_z).T - 273.15
print("31.进口蒸汽温度", format(Tl_i, ".2f"))
delta_pcd = 0.05 * pcd  # 低压缸排气压损
pl_z = pcd + delta_pcd
print("32.排汽压力", format(pl_z, ".4f"))
sl_i = IAPWS97(P=pl_i, h=hrh2_z).s
h_ilz = IAPWS97(P=pl_z, s=sl_i).h
hl_z = hrh2_z - eta_l_i * (hrh2_z - h_ilz)
xl_z = IAPWS97(P=pl_z, h=hl_z).x
print("33.排汽干度", format(xl_z, ".4f"))

Z = 7
print("34.回热级数", format(Z, ".0f"))
Zl = 4
print("35.低压给水加热器级数", format(Zl, ".0f"))
Zh = 2
print("36.高压给水加热器级数", format(Zh, ".0f"))

# 给水的焓升分配
h_s = IAPWS97(P=ps, x=0).h
hcd = IAPWS97(T=Tcd+273.15, x=0).h
delta_hfw_op = (h_s - hcd) / (Z + 1)
hfw_op = hcd + Z * delta_hfw_op
Tfw_op = IAPWS97(P=ps, h=hfw_op).T - 273.15
Tfw = 0.875 * Tfw_op    # 0.85-0.90
hfw = IAPWS97(P=ps, T=Tfw+273.15).h
delta_hfw = (hfw - hcd) / Z
print("37.第一次给水回热分配", format(delta_hfw, ".2f"))

# 第二次给水回热分配
pdea = 0.99 * ph_z
hdea_o = IAPWS97(P=pdea, x=0).h
delta_hfw_h = (hfw - hdea_o) / Zh
print("38.高压加热器给水焓升", format(delta_hfw_h, ".2f"))
delta_hfw_l = (hdea_o - hcd) / (Zl + 1)
print("39.除氧器及低加给水焓升", format(delta_hfw_l, ".2f"))

# 给水回热系统中的压力选择
pcwp = pdea * 3.1   # 3-3.2
hcwp = hcd  # 凝水泵出口给水比焓
delta_pcws = pcwp - pdea
delta_pfi = delta_pcws / (Zl + 1)   # 每级低压加热器及除氧器的平均压降

# 低压加热器给水参数

# j = 1
# plfwi_1 = pcwp  # 进口给水压力
# hlfwi_1 = hcwp  # 进口给水比焓
# Tlfwi_1 = IAPWS97(P=plfwi_1, h=hlfwi_1).T - 273.15  # 进口给水温度
# plfwo_1 = plfwi_1 - delta_pfi   # 出口给水压力
# hlfwo_1 = hlfwi_1 + delta_hfw_l  # 出口给水比焓
# Tlfwo_1 = IAPWS97(P=plfwo_1, h=hlfwo_1).T - 273.15  # 出口给水温度
# Tlrok_1 = Tlfwo_1 + theta_l_u  # 汽侧疏水温度
# hlrok_1 = IAPWS97(T=Tlrok_1+273.15, x=0).h  # 汽侧疏水比焓
# plrok_1 = IAPWS97(T=Tlrok_1+273.15, x=0).P  # 汽侧压力
# print("第{}级进口给水比焓".format(j), format(hlfwi_1, ".2f"))
# print("第{}级出口给水比焓".format(j), format(hlfwo_1, ".2f"))
# print("第{}级进口给水温度".format(j), format(Tlfwi_1, ".2f"))
# print("第{}级进口给水温度".format(j), format(Tlfwo_1, ".2f"))

j = 1
plfwi = np.zeros(5)
hlfwi = np.zeros(5)
Tlfwi = np.zeros(5)
plfwo = np.zeros(5)
hlfwo = np.zeros(5)
Tlfwo = np.zeros(5)
Tlrok = np.zeros(5)
hlrok = np.zeros(5)
plrok = np.zeros(5)
plfwi[j] = pcwp  # 进口给水压力
hlfwi[j] = hcwp  # 进口给水比焓
Tlfwi[j] = IAPWS97(P=plfwi[j], h=hlfwi[j]).T - 273.15  # 进口给水温度
plfwo[j] = plfwi[j] - delta_pfi   # 出口给水压力
hlfwo[j] = hlfwi[j] + delta_hfw_l  # 出口给水比焓
Tlfwo[j] = IAPWS97(P=plfwo[j], h=hlfwo[j]).T - 273.15  # 出口给水温度
Tlrok[j] = Tlfwo[j] + theta_l_u  # 汽侧疏水温度
hlrok[j] = IAPWS97(T=Tlrok[j]+273.15, x=0).h  # 汽侧疏水比焓
plrok[j] = IAPWS97(T=Tlrok[j]+273.15, x=0).P  # 汽侧压力
print("第{}级进口给水比焓".format(j), format(hlfwi[j], ".2f"))
print("第{}级出口给水比焓".format(j), format(hlfwo[j], ".2f"))
print("第{}级进口给水温度".format(j), format(Tlfwi[j], ".2f"))
print("第{}级进口给水温度".format(j), format(Tlfwo[j], ".2f"))

while j <= 3:
    j = j + 1
    plfwi[j] = plfwo[j-1]
    hlfwi[j] = hlfwo[j-1]
    Tlfwi[j] = Tlfwo[j-1]
    plfwo[j] = plfwi[j] - delta_pfi
    hlfwo[j] = hlfwi[j] + delta_hfw_l
    Tlfwo[j] = IAPWS97(P=plfwo[j], h=hlfwo[j]).T - 273.15  # 出口给水温度
    Tlrok[j] = Tlfwo[j] + theta_l_u  # 汽侧疏水温度
    hlrok[j] = IAPWS97(T=Tlrok[j] + 273.15, x=0).h  # 汽侧疏水比焓
    plrok[j] = IAPWS97(T=Tlrok[j] + 273.15, x=0).P  # 汽侧压力
    print("第{}级进口给水比焓".format(j), format(hlfwi[j], ".2f"))
    print("第{}级出口给水比焓".format(j), format(hlfwo[j], ".2f"))
    print("第{}级进口给水温度".format(j), format(Tlfwi[j], ".2f"))
    print("第{}级进口给水温度".format(j), format(Tlfwo[j], ".2f"))

# 除氧器
hdea_i = hlfwo[4]
print("41.进口给水比焓", format(hdea_i, ".2f"))
hdea_o = hdea_i + delta_hfw
print("42.出口给水比焓", format(hdea_o, ".2f"))
Tdea = IAPWS97(P=pdea, x=0).T - 273.15
print("43.出口给水温度", format(Tdea, ".2f"))
print("44.运行压力", format(pdea, ".2f"))
pfwpo = 1.2 * ps    # 1.15-1.25
hfwpo = hdea_o  # 给水泵出口流体比焓
pfwi = ps + 0.1     # GS二次侧进口给水压力

# 高压加热器给水参数
# 第六级高压给水加热器
phfwi_6 = pfwpo
hhfwi_6 = hfwpo
print("第6级进口给水比焓", format(hhfwi_6, ".2f"))
Thfwi_6 = IAPWS97(P=phfwi_6, h=hhfwi_6).T - 273.15
print("第6级进口给水温度", format(Thfwi_6, ".2f"))
phfwo_6 = phfwi_6 - (phfwi_6 - pfwi) / 2
hhfwo_6 = hhfwi_6 + delta_hfw_h
print("第6级出口给水比焓", format(hhfwo_6, ".2f"))
Thfwo_6 = IAPWS97(P=phfwo_6, h=hhfwo_6).T - 273.15
print("第6级出口给水温度", format(Thfwo_6, ".2f"))
Throk_6 = Thfwo_6 + theta_h_u  # 汽侧疏水温度
hhrok_6 = IAPWS97(T=Throk_6+273.15, x=0).h
phrok_6 = IAPWS97(T=Throk_6+273.15, x=0).P

# 第七级高压给水加热器
phfwi_7 = phfwo_6
hhfwi_7 = hhfwo_6
print("第7级进口给水比焓", format(hhfwi_7, ".2f"))
Thfwi_7 = Thfwo_6
print("第7级进口给水温度", format(Thfwi_7, ".2f"))
phfwo_7 = pfwi
hhfwo_7 = hhfwi_7 + delta_hfw_h
print("第7级出口给水比焓", format(hhfwo_7, ".2f"))
Thfwo_7 = IAPWS97(P=phfwo_7, h=hhfwo_7).T - 273.15
print("第7级出口给水温度", format(Thfwo_7, ".2f"))
Throk_7 = Thfwo_7 + theta_h_u
hhrok_7 = IAPWS97(T=Throk_7+273.15, x=0).h
phrok_7 = IAPWS97(T=Throk_7+273.15, x=0).P

# 高压缸抽汽
delta_pe_j = 0.04   # 0.03-0.05

# 第六级给水加热器抽汽参数
delta_pe_6 = delta_pe_j * phrok_6  # 压损
phes_6 = phrok_6 + delta_pe_6
print("第6级抽汽压力", format(phes_6, ".2f"))
hihes_6 = IAPWS97(P=phes_6, s=sh_i).h  # 抽汽理想比焓
hhes_6 = hhi - eta_m * eta_h_i * (hhi - hihes_6)  # 抽汽比焓
xhes_6 = IAPWS97(P=phes_6, h=hhes_6).x  # 抽汽干度
print("第6级抽汽干度", format(xhes_6, ".4f"))
# 第七级给水加热器抽汽参数
delta_pe_7 = delta_pe_j * phrok_7  # 压损
phes_7 = phrok_7 + delta_pe_7
print("第7级抽汽压力", format(phes_7, ".2f"))
hihes_7 = IAPWS97(P=phes_7, s=sh_i).h  # 抽汽理想比焓
hhes_7 = hhi - eta_m * eta_h_i * (hhi - hihes_7)  # 抽汽比焓
xhes_7 = IAPWS97(P=phes_7, h=hhes_7).x  # 抽汽干度
print("第7级抽汽干度", format(xhes_7, ".4f"))

# 低压缸抽汽
j = 1
delta_pe = np.zeros(5)
ples = np.zeros(5)
hiles = np.zeros(5)
hles = np.zeros(5)
xles = np.zeros(5)
while j <= 4:
    delta_pe[j] = delta_pe_j * plrok[j]
    ples[j] = plrok[j] + delta_pe[j]
    print("第{}级抽汽压力".format(j), format(ples[j], ".2f"))
    hiles[j] = IAPWS97(P=ples[j], s=sl_i).h
    hles[j] = hli - eta_m * eta_l_i * (hli - hiles[j])
    xles[j] = IAPWS97(P=ples[j], h=hles[j]).x
    print("第{}级抽汽干度".format(j), format(xles[j], ".4f"))
    j = j + 1

# 再热器抽汽
# 一级再热器抽汽参数
prh_1 = phes_7  # 加热蒸汽进口压力
xrh_1 = xhes_7  # 加热蒸汽进口干度
Trh_1 = IAPWS97(P=prh_1, x=xrh_1).T - 273.15  # 加热蒸汽进口温度
hrh_1 = IAPWS97(P=prh_1, x=xrh_1).h  # 加热蒸汽进口比焓
hzs_1 = IAPWS97(P=prh_1, x=0).h  # 再热器疏水比焓
# print(prh_1,xrh_1,Trh_1,hrh_1,hzs_1)
# 二级再热器抽汽参数
prh_2 = ph_i
xrh_2 = xh_i
Trh_2 = IAPWS97(P=prh_2, x=xrh_2).T - 273.15  # 加热蒸汽进口温度
hrh_2 = IAPWS97(P=prh_2, x=xrh_2).h  # 加热蒸汽进口比焓
hzs_2 = IAPWS97(P=prh_2, x=0).h  # 再热器疏水比焓
# print(prh_2,xrh_2,Trh_2,hrh_2,hzs_2)

# 参数校核
if 20 <= delta_Tm <= 33:
    print("一、二次侧对数平均温差校核通过")
else:
    print("一、二次侧对数平均温差校核失败")

if 4.2 * 10 ** (-3) <= pcd <= 7.5 * 10 ** (-3):
    print("冷凝器压力校核通过")
else:
    print("冷凝器压力校核失败")

if xh_z >= 0.84 and xl_z >= 0.84:
    print("高、低压缸排汽湿度校核通过")
else:
    print("高、低压缸排汽湿度校核失败")


# 蒸汽发生器总蒸汽产量的计算
eta_e_NPP = 0.3     # 假定电厂效率
Gcd = 2000     # 假定冷凝器凝水量
Ha = hhi - IAPWS97(P=1.05*pcd, x=0).h  # 给水泵汽轮机中蒸汽的绝热焓降
for i in range(0, 1000):
    QR = Ne / eta_e_NPP  # 反应堆功率
    Ds = 1000 * QR * eta_1 / ((hfh - h_s) + (1 + xi_d) * (h_s - hfw))  # GS的蒸汽产量
    Gfw = (1 + xi_d) * Ds   # 给水流量
    for j in range(0, 1000):
        # 给水泵的扬程
        Hfwp = pfwpo - pdea
        # 给水泵中给水的密度,定为给水泵密度进出口平均值
        rho_fw = (IAPWS97(P=pdea, x=0).rho + IAPWS97(P=pfwpo, x=0).rho) / 2
        # 给水泵有效输出功率
        Nfwpp = Gfw * Hfwp / rho_fw
        # 给水泵汽轮机耗气量
        Nfwp_t = Nfwpp / (eta_fwp_p * eta_fwp_ti * eta_fwp_tm * eta_fwp_tg)     # 给水泵汽轮机理论功率
        Gs_fwp = 1000 * Nfwp_t / Ha
        # 低压给水加热器抽汽量
        Gles_4 = Gcd * (hlfwo[4] - hlfwi[4]) / (eta_h * (hles[4] - hlrok[4]))
        Gles_3 = ((Gcd * (hlfwo[3] - hlfwi[3]) - eta_h * Gles_4 * (hlrok[4] - hlrok[3]))
                  / (eta_h * (hles[3] - hlrok[3])))
        Gles_2 = ((Gcd * (hlfwo[2] - hlfwi[2]) - eta_h * (Gles_3 + Gles_4) * (hlrok[3] - hlrok[2]))
                  / (eta_h * (hles[2] - hlrok[2])))
        Gles_1 = ((Gcd * (hlfwo[1] - hlfwi[1]) - eta_h * (Gles_2 + Gles_3 + Gles_4) * (hlrok[2] - hlrok[1]))
                  / (eta_h * (hles[1] - hlrok[1])))
        # 低压缸耗气量
        Gs_lp = (0.6 * 1000 * Ne / (eta_m * eta_ge) + Gles_4 * (hles[4] - hl_z) + Gles_3 * (hles[3] - hl_z)
                 + Gles_2 * (hles[2] - hl_z) + Gles_1 * (hles[1] - hl_z)) / (hli-hl_z)
        # 再热器加热蒸汽量
        Gs_rh_1 = Gs_lp * delta_hrh / (eta_h * (hrh_1 - hzs_1))
        Gs_rh_2 = Gs_lp * delta_hrh / (eta_h * (hrh_2 - hzs_2))
        # 高压给水加热器抽汽量
        Ghes_7 = ((Gfw * (hhfwo_7 - hhfwi_7) - eta_h * Gs_rh_2 * (hzs_2 - hhrok_7))
                  / (eta_h * (hhes_7 - hhrok_7)))
        Ghes_6 = ((Gfw * (hhfwo_6 - hhfwi_6) - eta_h * Gs_rh_1 * (hzs_1 - hhrok_6) - eta_h * (Gs_rh_2 + Ghes_7) * (hhrok_7 - hhrok_6))
                  / (eta_h * (hhes_6 - hhrok_6)))
        # 汽水分离器疏水流量
        Guw = Gs_lp * (xrh1_i - xsp_i) / xsp_i
        G_h1 = Gs_lp + Guw
        # 除氧器耗气量
        Gs_dea = (Gfw * hdea_o - Guw * huw - Gcd * hlfwo[4] - (Gs_rh_1 + Gs_rh_2 + Ghes_6 + Ghes_7) * hhrok_6) / hh_z
        # 高压缸出口排气总流量
        Gh = G_h1 + Gs_lp + Guw
        # 高压缸耗气量
        Gs_hp = (0.4 * 1000 * Ne / (eta_m * eta_ge) + Ghes_7 * (hhes_7 - hh_z) + Ghes_6 * (hhes_6 - hh_z)
                 + Gs_rh_1 * (hrh_1 - hh_z))/(hhi - hh_z)
        # 对假设冷凝水流量验证
        Ds = Gs_fwp + Gs_rh_2 + Gs_hp  # 对应的新蒸汽耗量
        Gfw1 = (1 + xi_d) * Ds  # 给水流量
        Gcd1 = Gfw1 - Gs_dea - Guw - (Ghes_6 + Ghes_7 + Gs_rh_1 + Gs_rh_2)
        if abs(Gcd1 - Gcd) / Gcd < 1e-2:
            break
        else:
            Gcd = Gcd1
            Gfw = Gfw1
    QR1 = (Ds * (hfh - hfw) + xi_d * Ds * (h_s - hfw)) / (1000 * eta_1)
    eta_e_NPP1 = Ne / QR1

    # 输出结果
    print('1.核电厂效率', format(eta_e_NPP, ".4f"))
    print('2.反应堆热功率', format(QR, ".2f"))
    print('3.蒸汽发生器总蒸汽产量', format(Ds, ".2f"))
    print('4.汽轮机高压缸耗汽量', format(Gs_hp, ".2f"))
    print('5.汽轮机低压缸耗汽量', format(Gs_lp, ".2f"))
    print('6.第一级再热器耗汽量', format(Gs_rh_1, ".2f"))
    print('7.第二级再热器耗汽量', format(Gs_rh_2, ".2f"))
    print('8.除氧器耗汽量', format(Gs_dea, ".2f"))
    print('9.给水泵汽轮机耗汽量', format(Gs_fwp, ".2f"))
    print('10.给水泵给水量', format(Gfw, ".2f"))
    print('11.给水泵扬程', format(Hfwp, ".2f"))
    print('12.高压缸抽汽量')
    print(' 12.1.第七级抽汽量', format(Ghes_7, ".2f"))
    print(' 12.2.第六级抽汽量', format(Ghes_6, ".2f"))
    print('13.低压缸抽汽量')
    print(' 13.1.第四级抽汽量', format(Gles_4, ".2f"))
    print(' 13.2.第三级抽汽量', format(Gles_3, ".2f"))
    print(' 13.3.第二级抽汽量', format(Gles_2, ".2f"))
    print(' 13.4.第一级抽汽量', format(Gles_1, ".2f"))
    print('14.凝结水量', format(Gcd, ".2f"))
    print('15.汽水分离器疏水量', format(Guw, ".2f"))

    if abs(eta_e_NPP - eta_e_NPP1) / eta_e_NPP < 1e-3:
        break
    else:
        eta_e_NPP = eta_e_NPP1
        print("==========继续迭代==========")
print("==========迭代完成==========")
