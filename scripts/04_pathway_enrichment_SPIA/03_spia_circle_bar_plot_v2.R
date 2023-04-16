df <- unique(df2)

df <- rbind(df %>% filter(tA > 0) %>% arrange(desc(tA)),
             df %>% filter(tA < 0) %>% arrange(tA))

df$pos1 <- 1:nrow(df)

df$tA <- ifelse(abs(df$tA) > 0.1, df$tA, sign(df$tA)*0.1)

df %>%
  arrange(tA) %>%
  filter(tA > 0) %>%
  mutate(Name=factor(Name, Name)) %>%
  ggplot( aes(x=Name, y=tA) ) +
  geom_bar(stat="identity", fill = "#cb4154") +
  geom_text(hjust = 1, size = 10, fontface = "bold", 
            aes( y = 0, label = paste(Name," "))) +
  # theme_ipsum() +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position="none",
    axis.text = element_blank()
  ) +
  xlab("") +
  ylab("") +
  coord_polar(theta = "y", clip = "off") +
  ylim(0,ceiling(max(abs(range(df$tA)))) + 3) +
  theme_minimal()+
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(colour = "black", size = 30))

wd <- df %>%
  arrange(tA) %>%
  filter(tA > 0) %>% 
  nrow()

ggsave("circleBarClust4IDC.png", width = wd + 5, height = wd, units = "in", bg = "white")

df %>%
  arrange(desc(tA)) %>%
  filter(tA < 0) %>%
  mutate(Name=factor(Name, Name)) %>%
  ggplot( aes(x=Name, y=-tA) ) +
  geom_bar(stat="identity", fill = "skyblue") +
  geom_text(hjust = 0, size = 10, fontface = "bold", 
            aes( y = 0, label = paste(" ", Name))) +
  # theme_ipsum() +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position="none",
    axis.text = element_blank()
  ) +
  xlab("") +
  ylab("") +
  coord_polar(theta = "y", clip = "off") +
  ylim(ceiling(max(abs(range(df$tA)))) + 3, 0) +
  theme_minimal()+
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(colour = "black", size = 30))

sd <- df %>%
  arrange(tA) %>%
  filter(tA < 0) %>% 
  nrow()

ggsave("circleBarClust4ILC.png", width = wd + 5, height = wd, units = "in", bg = "white")

png("../3.Result/spiaVisulBarCir2.png", width =25, height = 10, units = "in", res = 300)
ggarrange(p1, p2, ncol=2, nrow=1, common.legend = TRUE, legend="bottom", widths = c(11, 4))
dev.off()

df %>%
  # arrange(ifelse(tA > 0, desc(tA), tA)) %>%
  # filter(tA < 0) %>%
  mutate(Name=factor(Name, Name)) %>%
  ggplot( aes(x=-pos1, y=tA, fill = Status) ) +
  geom_bar(stat="identity", width = 0.95) +
  geom_text(hjust = ifelse(df$tA > 0, 1, 0), size = 5, fontface = "bold", 
            aes( y = 0, 
                 label = ifelse(tA > 0, paste(Name, " "), paste(" ",Name)))) +
  # theme_ipsum() +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position="none",
    axis.text = element_blank()
  ) +
  xlab("") +
  ylab("") +
  coord_polar(theta = "y", start = 18.65/4, clip = "on") +
  ylim(range(df$tA)) +
  expand_limits(x = 0)
ggsave("test.png", width = 15, height = 15, units = "in")
