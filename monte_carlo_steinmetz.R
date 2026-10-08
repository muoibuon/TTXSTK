# =====================================================================
# Chủ đề 10: Tích phân Monte Carlo trong Hình học
# Thí nghiệm: Ước lượng thể tích khối Steinmetz (giao của 2 hình trụ)
#   Hình trụ 1: x^2 + y^2 <= 1     Hình trụ 2: x^2 + z^2 <= 1
#   Thể tích chính xác: V = 16/3 ≈ 5.3333
# Công cụ: Phân bố đều U(-1, 1) + Luật Yếu Số Lớn (WLLN)
# =====================================================================

set.seed(2026)                       # cố định kết quả để tái lập được

V_true <- 16 / 3                     # đáp số lý thuyết để đối chiếu
V_box  <- 2^3                        # thể tích hộp bao [-1,1]^3 = 8

# ---------------------------------------------------------------------
# Bước 1: Gieo N điểm ngẫu nhiên đều trong hộp [-1,1]^3
#   Mỗi tọa độ x, y, z ~ U(-1, 1) độc lập
# ---------------------------------------------------------------------
N <- 100000
x <- runif(N, min = -1, max = 1)
y <- runif(N, min = -1, max = 1)
z <- runif(N, min = -1, max = 1)

# ---------------------------------------------------------------------
# Bước 2: Biến chỉ thị I_i = 1 nếu điểm nằm trong khối, 0 nếu không
#   I_i ~ Bernoulli(p), với p = V / V_box
# ---------------------------------------------------------------------
inside <- (x^2 + y^2 <= 1) & (x^2 + z^2 <= 1)

# ---------------------------------------------------------------------
# Bước 3: Ước lượng  V_hat = V_box * (1/N) * sum(I_i)
#   Theo WLLN: trung bình mẫu -> p (hội tụ theo xác suất) khi N -> ∞
# ---------------------------------------------------------------------
p_hat <- mean(inside)
V_hat <- V_box * p_hat

# Sai số chuẩn và khoảng tin cậy 95%
se    <- V_box * sqrt(p_hat * (1 - p_hat) / N)
ci    <- V_hat + c(-1, 1) * 1.96 * se

cat("===== KẾT QUẢ (N =", format(N, big.mark = ",", scientific = FALSE), "điểm) =====\n")
cat(sprintf("Tỷ lệ điểm rơi vào khối p_hat : %.5f\n", p_hat))
cat(sprintf("Thể tích ước lượng V_hat      : %.5f\n", V_hat))
cat(sprintf("Thể tích chính xác 16/3       : %.5f\n", V_true))
cat(sprintf("Sai số tuyệt đối              : %.5f\n", abs(V_hat - V_true)))
cat(sprintf("Khoảng tin cậy 95%%            : [%.4f, %.4f]\n", ci[1], ci[2]))

# ---------------------------------------------------------------------
# Bước 4: Minh họa WLLN — ước lượng hội tụ khi N tăng
# ---------------------------------------------------------------------
running_est <- V_box * cumsum(inside) / seq_len(N)

par(mfrow = c(1, 2))

plot(seq_len(N), running_est, type = "l", log = "x", col = "steelblue",
     ylim = c(4.5, 6.5), xlab = "Số điểm N (thang log)",
     ylab = "Thể tích ước lượng",
     main = "Hội tụ của ước lượng (WLLN)")
abline(h = V_true, col = "red", lwd = 2, lty = 2)
legend("topright", legend = c("Ước lượng", "Giá trị đúng 16/3"),
       col = c("steelblue", "red"), lty = c(1, 2), bty = "n", cex = 0.8)

# ---------------------------------------------------------------------
# Bước 5: Kiểm tra định lượng — xác suất sai lệch > epsilon giảm về 0
#   Lặp lại thí nghiệm R lần với mỗi cỡ mẫu n
# ---------------------------------------------------------------------
mc_volume <- function(n) {
  x <- runif(n, -1, 1); y <- runif(n, -1, 1); z <- runif(n, -1, 1)
  V_box * mean((x^2 + y^2 <= 1) & (x^2 + z^2 <= 1))
}

n_vals <- c(100, 1000, 10000, 100000)
R      <- 200
eps    <- 0.1

res <- t(sapply(n_vals, function(n) {
  est <- replicate(R, mc_volume(n))
  c(n          = n,
    trung_binh = mean(est),
    do_lech_sd = sd(est),
    P_sai_lech = mean(abs(est - V_true) > eps))
}))

cat("\n===== WLLN: P(|V_hat - V| >", eps, ") khi n tăng =====\n")
print(round(res, 4))

plot(res[, "n"], res[, "do_lech_sd"], log = "xy", type = "b", pch = 19,
     col = "darkgreen", xlab = "Cỡ mẫu n", ylab = "Độ lệch chuẩn",
     main = "Sai số giảm theo 1/sqrt(n)")
lines(n_vals, res[1, "do_lech_sd"] * sqrt(n_vals[1] / n_vals),
      col = "red", lty = 2)
legend("topright", legend = c("Thực nghiệm", "Lý thuyết ~ 1/sqrt(n)"),
       col = c("darkgreen", "red"), lty = c(1, 2), pch = c(19, NA),
       bty = "n", cex = 0.8)

par(mfrow = c(1, 1))
