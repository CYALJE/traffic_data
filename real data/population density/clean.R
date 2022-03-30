dt2 = read.csv("Documents/real data/population density/pd2.csv", header=T)
dt3 = read.csv("Documents/real data/population density/pd3.csv", header=T)
dt4 = read.csv("Documents/real data/population density/pd4.csv", header=T)
dt5 = read.csv("Documents/real data/population density/pd5.csv", header=T)
dt6 = read.csv("Documents/real data/population density/pd6.csv", header=T)
dt7 = read.csv("Documents/real data/population density/pd7.csv", header=T)
dt8 = read.csv("Documents/real data/population density/pd8.csv", header=T)
dt9 = read.csv("Documents/real data/population density/pd9.csv", header=T)
dt10 = read.csv("Documents/real data/population density/pd10.csv", header=T)


total = rbind(dt2, dt3, dt4, dt5, dt6,
              dt7, dt8, dt9, dt10)
head(total)
con = c("新北市永和區", "新北市蘆洲區", "臺北市大安區", "高雄市新興區", "新北市板橋區",
        "新北市三重區", "新北市新莊區", "臺北市大同區", "臺北市松山區", "臺中市北區",
        "高雄市苓雅區", "新北市中和區", "臺中市中區", "臺北市萬華區", "臺北市中正區",
        "臺中市西區", "高雄市旗津區", "臺北市信義區", "臺中市南區", "高雄市三民區", 
        "高雄市鹽埕區")
select = total[(total[, 2] %in% con) == TRUE, ]
View(select)
write.csv(select, "Documents/real data/population density/clean.csv", fileEncoding="big-5")


site.order=c("高雄市苓雅區",
             "高雄市旗津區",
             "高雄市三民區",
             "高雄市新興區",
             "高雄市鹽埕區",
             "臺北市大安區",
             "臺北市大同區",
             "臺北市松山區",
             "臺北市萬華區",
             "臺北市信義區",
             "臺北市中正區",
             "臺中市北區",
             "臺中市南區",
             "臺中市西區",
             "臺中市中區",
             "新北市板橋區",
             "新北市蘆洲區",
             "新北市三重區",
             "新北市新莊區",
             "新北市永和區",
             "新北市中和區")

clean = read.csv("Documents/real data/population density/clean.csv", header=T, fileEncoding="big-5") 
head(clean)

clean.copy = clean
library(tidyverse)
tmp = clean.copy %>% arrange(statistic_yyy, site_id)
tmp = as.data.frame(tmp)
write.csv(tmp, "test.csv", fileEncoding="big-5")

