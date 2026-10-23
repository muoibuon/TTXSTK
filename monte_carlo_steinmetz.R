# =====================================================================
# Chủ đề 10: Tích phân Monte Carlo trong Hình học
# Thí nghiệm: Ước lượng thể tích khối Steinmetz (giao của 2 hình trụ)
#   Hình trụ 1: x^2 + y^2 <= 1     Hình trụ 2: x^2 + z^2 <= 1
#   Thể tích chính xác: V = 16/3 ≈ 5.3333
# Công cụ: Phân bố đều U(-1, 1) + Luật Yếu Số Lớn (WLLN)
#
# File chia 3 phần, khớp 3 mục trên trang web:
#   PHẦN 01 - Gieo điểm vào hộp [-1,1]^3, ước lượng V_hat, khoảng tin cậy
#   PHẦN 02 - Vì sao đáp số là 16/3: lát cắt hình vuông và tích phân
#   PHẦN 03 - Luật Yếu Số Lớn: lặp lại thí nghiệm, Chebyshev, 1/sqrt(n)
# =====================================================================

set.seed(2026)                       # cố định kết quả để tái lập được

V_true <- 16 / 3                     # đáp số lý thuyết để đối chiếu
V_box  <- 2^3                        # thể tích hộp bao [-1,1]^3 = 8
p_true <- V_true / V_box             # xác suất trúng p = 2/3

# Điểm (x, y, z) nằm trong khối khi nằm trong cả hai hình trụ
in_khoi <- function(x, y, z) (x^2 + y^2 <= 1) & (x^2 + z^2 <= 1)

par(mfrow = c(2, 2))                 # 4 đồ thị: 01, 02 (hai hình), 03

# #####################################################################
# PHẦN 01 - GIEO ĐIỂM VÀO HỘP [-1,1]^3
# #####################################################################

# ---------------------------------------------------------------------
# Bước 1: Gieo N điểm ngẫu nhiên đều trong hộp [-1,1]^3
#   Mỗi tọa độ x, y, z ~ U(-1, 1) độc lập, mật độ f = 1/2 trên [-1, 1]
#   => mật độ trong hộp là 1/8, nên P(điểm rơi vào khối) = V / 8
# ---------------------------------------------------------------------
N <- 100000
x <- runif(N, min = -1, max = 1)
y <- runif(N, min = -1, max = 1)
z <- runif(N, min = -1, max = 1)

# ---------------------------------------------------------------------
# Bước 2: Biến chỉ thị I_i = 1 nếu điểm nằm trong khối, 0 nếu không
#   I_i ~ Bernoulli(p), với p = V / V_box
# ---------------------------------------------------------------------
inside <- in_khoi(x, y, z)

# ---------------------------------------------------------------------
# Bước 3: Ước lượng  V_hat = V_box * (1/N) * sum(I_i)
#   Theo WLLN: trung bình mẫu -> p (hội tụ theo xác suất) khi N -> ∞
# ---------------------------------------------------------------------
p_hat <- mean(inside)
V_hat <- V_box * p_hat

# Sai số chuẩn và khoảng tin cậy 95%
se    <- V_box * sqrt(p_hat * (1 - p_hat) / N)
ci    <- V_hat + c(-1, 1) * 1.96 * se

cat("===== PHẦN 01: KẾT QUẢ (N =", format(N, big.mark = ",", scientific = FALSE), "điểm) =====\n")
cat(sprintf("Tỷ lệ điểm rơi vào khối p_hat : %.5f\n", p_hat))
cat(sprintf("Thể tích ước lượng V_hat      : %.5f\n", V_hat))
cat(sprintf("Thể tích chính xác 16/3       : %.5f\n", V_true))
cat(sprintf("Sai số tuyệt đối              : %.5f\n", abs(V_hat - V_true)))
cat(sprintf("Sai số chuẩn                  : %.5f\n", se))
cat(sprintf("Khoảng tin cậy 95%%            : [%.4f, %.4f]\n", ci[1], ci[2]))

# Ví dụ kiểm tra từng điểm (giống bảng "Ví dụ" dưới mô phỏng 01 trên web)
vd <- data.frame(x = c(0.5, 0.9, 0.6, -0.8, 0.2, 0.95, 0.6, 0),
                 y = c(0.5, 0.9, 0.3, 0.7, -0.6, 0.1, 0.8, 1),
                 z = c(0.5, 0, 0.9, -0.7, 0.7, 0.1, 0, 1))
vd$ong1     <- vd$x^2 + vd$y^2
vd$ong2     <- vd$x^2 + vd$z^2
vd$trong_khoi <- in_khoi(vd$x, vd$y, vd$z)
cat("\nVí dụ kiểm tra điểm (trong khối khi cả ong1 và ong2 đều <= 1):\n")
print(vd, row.names = FALSE)

# ---------------------------------------------------------------------
# Bước 4: Minh họa WLLN - ước lượng hội tụ khi N tăng
# ---------------------------------------------------------------------
running_est <- V_box * cumsum(inside) / seq_len(N)
n_seq <- seq_len(N)
band  <- 1.96 * V_box * sqrt(p_true * (1 - p_true) / n_seq)

plot(n_seq, running_est, type = "l", log = "x", col = "steelblue",
     ylim = c(4.5, 6.5), xlab = "Số điểm N (thang log)",
     ylab = "Thể tích ước lượng",
     main = "01. Hội tụ của ước lượng (WLLN)")
lines(n_seq, V_true + band, col = "grey60", lty = 3)
lines(n_seq, V_true - band, col = "grey60", lty = 3)
abline(h = V_true, col = "red", lwd = 2, lty = 2)
legend("topright", legend = c("Ước lượng", "Giá trị đúng 16/3", "Vùng 95% lý thuyết"),
       col = c("steelblue", "red", "grey60"), lty = c(1, 2, 3), bty = "n", cex = 0.8)

# #####################################################################
# PHẦN 02 - VÌ SAO ĐÁP SỐ LÀ 16/3 (kiểm tra bằng giải tích)
#   Không dùng số ngẫu nhiên, nên không làm đổi kết quả phần 01 và 03.
# #####################################################################

# ---------------------------------------------------------------------
# Bước 5: Lát cắt tại x cố định
#   Hình trụ 1: y^2 <= 1 - x^2  =>  -h <= y <= h, với h = sqrt(1 - x^2)
#   Hình trụ 2: z^2 <= 1 - x^2  =>  -h <= z <= h
#   => lát cắt là hình vuông cạnh 2h, diện tích A(x) = 4(1 - x^2)
# ---------------------------------------------------------------------
A <- function(x) 4 * (1 - x^2)                 # diện tích lát cắt
G <- function(x) 4 * x - 4 * x^3 / 3           # nguyên hàm của A(x)
V_tichluy <- function(x) G(x) - G(-1)          # thể tích từ -1 tới x

x_vd <- c(-1, -0.5, 0, 0.4, 0.6, 0.8, 1)
lat_cat <- data.frame(x         = x_vd,
                      mot_tru_x2 = 1 - x_vd^2,
                      h         = sqrt(1 - x_vd^2),
                      canh_2h   = 2 * sqrt(1 - x_vd^2),
                      A_x       = A(x_vd),
                      tich_luy  = V_tichluy(x_vd))

cat("\n===== PHẦN 02: LÁT CẮT VÀ TÍCH PHÂN =====\n")
print(round(lat_cat, 4), row.names = FALSE)

# ---------------------------------------------------------------------
# Bước 6: Cộng các lát = tích phân
#   V = ∫_{-1}^{1} 4(1 - x^2) dx = G(1) - G(-1) = 8/3 - (-8/3) = 16/3
# ---------------------------------------------------------------------
tp <- integrate(A, -1, 1)
cat(sprintf("Tích phân bằng integrate()    : %.6f\n", tp$value))
cat(sprintf("Tính tay G(1) - G(-1)         : %.6f\n", G(1) - G(-1)))
cat(sprintf("So với V_hat ở phần 01        : %.5f (lệch %.5f)\n",
            V_hat, abs(V_hat - tp$value)))

# Hình lát cắt tại x = 0.4 (giống ô bên trái mục 02 trên web)
x0 <- 0.4; h0 <- sqrt(1 - x0^2)
gan <- abs(x - x0) < 0.03                      # các điểm đã gieo gần lát cắt
plot(y[gan], z[gan], pch = 16, cex = 0.4, asp = 1,
     col = ifelse(inside[gan], "seagreen3", "grey70"),
     xlim = c(-1, 1), ylim = c(-1, 1), xlab = "y", ylab = "z",
     main = sprintf("02. Lát cắt x = %.1f, cạnh %.3f", x0, 2 * h0))
rect(-h0, -h0, h0, h0, border = "seagreen4", lwd = 2)

# Đường cong A(x), tô phần thể tích tích lũy từ -1 tới x0
curve(A, -1, 1, lwd = 2, col = "seagreen4", xlab = "x", ylab = "A(x)",
      main = "02. Diện tích lát A(x) = 4(1 - x^2)")
xs <- seq(-1, x0, length.out = 100)
polygon(c(-1, xs, x0), c(0, A(xs), 0), col = adjustcolor("seagreen3", 0.3), border = NA)
segments(x0, 0, x0, A(x0), lty = 2)
legend("topright", bty = "n", cex = 0.8,
       legend = sprintf("Tích lũy tới x = %.1f: %.4f", x0, V_tichluy(x0)))

# #####################################################################
# PHẦN 03 - LUẬT YẾU SỐ LỚN: LẶP LẠI THÍ NGHIỆM
#   P(|V_hat - V| > eps) -> 0 khi n -> ∞
# #####################################################################

# ---------------------------------------------------------------------
# Bước 7: Lặp lại thí nghiệm R lần với mỗi cỡ mẫu n
# ---------------------------------------------------------------------
mc_volume <- function(n) {
  x <- runif(n, -1, 1); y <- runif(n, -1, 1); z <- runif(n, -1, 1)
  V_box * mean(in_khoi(x, y, z))
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

# ---------------------------------------------------------------------
# Bước 8: So với lý thuyết
#   SD lý thuyết  = 8 * sqrt(p(1-p)/n) ≈ 3.7712 / sqrt(n)
#   Cận Chebyshev = Var(V_hat) / eps^2 = 64 p(1-p) / (n eps^2), tối đa 1
# ---------------------------------------------------------------------
var_1diem <- V_box^2 * p_true * (1 - p_true)   # 64 * 2/9 ≈ 14.22
sd_ly_thuyet <- sqrt(var_1diem / n_vals)
can_cheb     <- pmin(1, var_1diem / (n_vals * eps^2))
res <- cbind(res, sd_ly_thuyet = sd_ly_thuyet, can_Chebyshev = can_cheb)

cat("\n===== PHẦN 03: WLLN, P(|V_hat - V| >", eps, ") khi n tăng =====\n")
print(round(res, 4))

# Tính ngược: cần bao nhiêu điểm
n_can_sd   <- ceiling(var_1diem / 0.01^2)             # để SD <= 0.01
n_can_cheb <- ceiling(var_1diem / (0.05 * eps^2))      # để cận Chebyshev <= 0.05
cat(sprintf("Cần n >= %s điểm để sai số chuẩn <= 0.01\n", format(n_can_sd, big.mark = ",")))
cat(sprintf("Cần n >= %s điểm để cận Chebyshev <= 0.05 (eps = %.1f)\n",
            format(n_can_cheb, big.mark = ","), eps))

plot(res[, "n"], res[, "do_lech_sd"], log = "xy", type = "b", pch = 19,
     col = "darkgreen", xlab = "Cỡ mẫu n", ylab = "Độ lệch chuẩn",
     main = "03. Sai số giảm theo 1/sqrt(n)")
lines(n_vals, sd_ly_thuyet, col = "red", lty = 2)
legend("topright", legend = c("Thực nghiệm", "Lý thuyết 3.77/sqrt(n)"),
       col = c("darkgreen", "red"), lty = c(1, 2), pch = c(19, NA),
       bty = "n", cex = 0.8)

par(mfrow = c(1, 1))
