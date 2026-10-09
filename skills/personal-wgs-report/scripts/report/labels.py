"""AADR 群体标签 → 中文说明（博客用）。先查手工词典，查不到就按“地区_遗址_时期”规则拼。"""
import re

MANUAL = {
    "China_Shaanxi_Wuzhuangguoliang_LN": "陕西五庄果墚 · 新石器晚期",
    "China_Henan_Haojiatai_IA": "河南郝家台 · 铁器时代",
    "China_Henan_Haojiatai_LN": "河南郝家台 · 龙山文化",
    "China_Shandong_Dinggong_LN": "山东丁公 · 龙山文化",
    "China_Shandong_Dinggong_LateAntiquity": "山东丁公 · 汉唐之间",
    "China_Weifang_XinZhi_Zhou": "山东潍坊 · 周代",
    "China_Henan_Pingliangtaisite_LN": "河南平粮台 · 龙山文化",
    "China_Zibo_TongLin_Jin": "山东淄博桐林 · 东周",
    "China_Zibo_YiXi": "山东淄博 · 东周",
    "China_Henan_Xisima_LShang": "河南西司马 · 晚商",
    "China_Henan_Xiaowusite_MN": "河南晓坞 · 仰韶文化",
    "China_Henan_Wanggousite_MN": "河南汪沟 · 仰韶文化",
    "Henan_Wadian": "河南瓦店 · 龙山文化",
    "China_Henan_Wadian_LN": "河南瓦店 · 龙山文化",
    "China_Henan_Wadiansite_LN": "河南瓦店 · 龙山文化",
    "Henan_Wangchenggang": "河南王城岗 · 龙山–二里头",
    "Henan_Yuzhuang": "河南禹州 · 龙山文化",
    "China_Henan_Jiaozuoniecunsite_LBA_IA": "河南焦作聂村 · 两周之际",
    "China_Shanxi_Shengedaliang_LN": "山西神圪垯梁 · 新石器晚期",
    "China_InnerMongolia_Miaozigousite_MN": "内蒙古庙子沟 · 仰韶晚期",
    "China_InnerMongolia_Erdaojingzi_LN": "内蒙古二道井子 · 夏家店下层",
    "China_Shandong_Chengziya_LN": "山东城子崖",
    "China_Shandong_Chengziya_Longshan": "山东城子崖 · 龙山文化",
    "China_Shandong_Chengziya_Yueshi": "山东城子崖 · 岳石文化",
    "China_Baligang_LN_Longshan": "河南八里岗 · 龙山文化",
    "China_Baligang_LN_Shijiahe": "河南八里岗 · 石家河文化",
    "China_Baligang_LN_Qujialing": "河南八里岗 · 屈家岭文化",
    "China_Baligang_N_Yangshao": "河南八里岗 · 仰韶文化",
    "China_Baligang_BA_EasternZhou": "河南八里岗 · 东周",
    "China_Jinan_LiuJiaZhuang_Shang": "山东济南刘家庄 · 商",
    "China_Qingdao_BeiQian_Dawenkou": "山东青岛北阡 · 大汶口文化",
    "China_Shandong_Bianbian_EN": "山东扁扁洞 · 新石器早期",
    "China_Shandong_BoshanMountain_EN": "山东博山 · 新石器早期",
    "China_TianyuanCave_UP": "北京田园洞人 · 4 万年前",
    "China_AmurRiverBasin_UP": "黑龙江流域 · 旧石器晚期",
    "China_AmurRiverBasin_LatePaleolithic": "黑龙江流域 · 旧石器晚期",
    "China_AmurRiverBasin_Mesolithic": "黑龙江流域 · 中石器",
    "China_AmurRiverBasin_N": "黑龙江流域 · 新石器",
    "China_Fujian_Tanshishan_LN": "福建昙石山 · 新石器晚期",
    "China_Fujian_Xitoucun_LN": "福建溪头村 · 新石器晚期",
    "China_Fujian_Qihedong_EN": "福建奇和洞 · 新石器早期",
    "China_Fujian_QiheCave_Epipaleolithic": "福建奇和洞 · 旧石器末期",
    "China_Guangxi_LonglinCave_Epipaleolithic": "广西隆林 · 旧石器末期",
    "China_Guangxi_BaojianshanCaveA_N": "广西宝剑山 · 新石器",
    "China_Guangxi_Dushan_N": "广西独山 · 新石器",
    "China_Penghu_Suogang_LN": "澎湖锁港 · 新石器晚期",
    "Taiwan_IA": "台湾汉本 · 铁器时代",
    "Taiwan_Gongguan": "台湾公馆",
    "China_Qinghai_Lajiasite_LN": "青海喇家 · 齐家文化",
    "China_Qinghai_Zongri": "青海宗日 · 新石器",
    "China_LN_Xiaoheyan": "内蒙古 · 小河沿文化",
    "China_Sitaimengguying_EN": "河北四台 · 新石器早期",
    "China_Xinjiang_Xiaohe_BA": "新疆小河 · 青铜时代",
    "Japan_Honshu_EarlyJomon": "日本绳文 · 早期",
    "Russia_PrimorskyKrai_AmurRiver_N": "俄罗斯滨海 · 新石器（鬼门洞）",
    "Mongolia_N": "蒙古 · 新石器",
    "Mongolia_East_N": "蒙古东部 · 新石器",
    "CHB": "北方汉族 · 北京（1000G）", "CHS": "南方汉族（1000G）", "CDX": "傣族 · 西双版纳（1000G）",
    "KHV": "越南京族（1000G）", "JPT": "日本 · 东京（1000G）", "Han": "汉族（HGDP）", "Japanese": "日本人（HGDP）",
    "Miao": "苗族", "She": "畲族", "Tujia": "土家族", "Yi": "彝族", "Naxi": "纳西族", "Dai": "傣族", "China_Lahu": "拉祜族",
    "Daur": "达斡尔族", "Hezhen": "赫哲族", "Oroqen": "鄂伦春族", "Xibo": "锡伯族", "Mongola": "蒙古族", "Tu": "土族",
    "Ami": "阿美人（台湾）", "Atayal": "泰雅人（台湾）", "Kinh_Vietnamese": "越南京族", "Thai": "泰人", "Uyghur": "维吾尔族",
    "Ulchi": "乌尔奇人", "Cambodian": "柬埔寨人", "Burmese": "缅甸人",
}
REGION = [("InnerMongolia", "内蒙古"), ("Henan", "河南"), ("Shandong", "山东"), ("Shaanxi", "陕西"), ("Shanxi", "山西"),
          ("Qinghai", "青海"), ("Xinjiang", "新疆"), ("Tibet", "西藏"), ("Guangxi", "广西"), ("Fujian", "福建"),
          ("Yunnan", "云南"), ("Sichuan", "四川"), ("Guizhou", "贵州"), ("Liaoning", "辽宁"), ("Jinan", "山东济南"),
          ("Zibo", "山东淄博"), ("Weifang", "山东潍坊"), ("Qingdao", "山东青岛"), ("Jining", "山东济宁"),
          ("Dongying", "山东东营"), ("Nagqu", "西藏那曲"), ("AmurRiverBasin", "黑龙江流域"), ("Mongolia", "蒙古"),
          ("Japan", "日本"), ("Taiwan", "台湾"), ("Vietnam", "越南"), ("Thailand", "泰国"), ("Laos", "老挝"),
          ("Russia", "俄罗斯"), ("Nepal", "尼泊尔"), ("Kazakhstan", "哈萨克斯坦"), ("Kyrgyzstan", "吉尔吉斯斯坦"),
          ("NorthernChina", "华北"), ("Penghu", "澎湖")]
PERIOD = [("LatePaleolithic", "旧石器晚期"), ("Epipaleolithic", "旧石器末期"), ("Mesolithic", "中石器"), ("UP", "旧石器晚期"),
          ("EN", "新石器早期"), ("MN", "新石器中期"), ("LN", "新石器晚期"), ("EBA", "青铜早期"), ("MBA", "青铜中期"),
          ("MLBA", "青铜中晚期"), ("LBA", "青铜晚期"), ("BA", "青铜时代"), ("EIA", "铁器早期"), ("IA", "铁器时代"),
          ("Historical", "历史时期"), ("Historic", "历史时期"), ("LateMedieval", "中世纪晚期"), ("EarlyMedieval", "中世纪早期"),
          ("Medieval", "中世纪"), ("XiongnuPeriod", "匈奴时期"), ("Xianbei", "鲜卑"), ("MingDynasty", "明代"),
          ("LateAntiquity", "汉唐之间"), ("Antiquity", "古典时期"), ("Shang", "商"), ("Zhou", "周"), ("Longshan", "龙山文化"),
          ("Yangshao", "仰韶文化"), ("Dawenkou", "大汶口文化"), ("Jomon", "绳文"), ("Yayoi", "弥生"), ("KofunPeriod", "古坟时代"),
          ("N", "新石器")]


def cn(group):
    if group in MANUAL:
        return MANUAL[group]
    g = re.sub(r"-\d+$", "", group)
    if g in MANUAL:
        return MANUAL[g]
    parts = g.split("_")
    if parts and parts[0] == "China":
        parts = parts[1:]
    reg = ""
    if parts:
        for k, v in REGION:
            if parts[0] == k:
                reg = v
                parts = parts[1:]
                break
    per = []
    rest = []
    for p in parts:
        hit = next((v for k, v in PERIOD if p == k), None)
        (per if hit else rest).append(hit or p)
    site = " ".join(re.sub(r"(site|Site|Cave|cave)$", "", r) for r in rest[:2])
    out = (reg + " " + site).strip()
    if per:
        out += " · " + "/".join(dict.fromkeys(per))
    return out or group
