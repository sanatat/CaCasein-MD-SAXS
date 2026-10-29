# نسخه ۴ — راهنمای فارسی

این نسخه یک **آزمون انتقادی روش** است؛ شبیه‌سازی واقعی MARTINI 3 نیست و ادعای اتصال واقعی کلسیم یا صحت تجربی SAXS ندارد.

## نتایج اصلی
- ۱۵ اجرای ۴۰۰۰ گامی: پنج حالت × سه تکرار مستقل.
- ۶ اجرای ۱۲۰۰۰ گامی برای ارزیابی اثر زمان و موقعیت اولیه یون.
- وقتی ۸ یون از سراسر جعبه به‌طور یکنواخت شروع شدند، میانگین تماس در نیمه دوم شبیه‌سازی کوتاه صفر بود؛ در شروع نزدیک به فسفوسرین حدود ۳٫۳۰ یون تماس داشتند.
- در آزمایش ۱۲۰۰۰ گامی نیز میانگین تماس در ثلث آخر برای آغاز دور از پروتئین **۰** و برای شروع نزدیک **۳٫۱۲** بود. این به معنای «عدم اتصال کلسیم واقعی» نیست؛ نشان می‌دهد نتایج به آغاز شبیه‌سازی وابسته‌اند و تعادل اثبات نشده است.

## راه‌اندازی Linux Mint

```bash
cd ~/Downloads
unzip CaCasein-MD-SAXS-v4-complete.zip
cd CaCasein-MD-SAXS
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python -m pytest -q tests v2_sequence_calcium/tests v3_explicit_ca/tests v4_validation/tests

cd v4_validation
python casein_v4.py --steps 4000 --replicas 3 --stride 100 --out new_results
python long_time_check.py --steps 12000 --replicas 3 --stride 120 --out new_long_results
```

برای سیستم تو (CPU بدون GPU) تنظیم پیش‌فرض مناسب است؛ در صورت کندی، `--steps 2000` و `--replicas 2` فقط برای تست سریع استفاده کن؛ نتایج گزارش اصلی با تنظیم‌های ذکرشده تولید شده‌اند.

## فایل‌ها
- `v4_validation/results/summary.csv`: آمار مستقل هر وضعیت؛
- `v4_validation/results/replica_summary.csv`: آمار هر تکرار؛
- `v4_validation/results/01_synthetic_saxs.png`: SAXS محاسبه‌شده در واحد فرضی؛
- `v4_validation/results/02_binding_control.png`: مقایسه تماس و پل‌های کلسیم؛
- `v4_validation/long_time_results/05_long_time_sampling.png`: اثر زمان و شروع تصادفی/نزدیک؛
- `v4_validation/RESULTS_V4_EN.md`: گزارش علمی انگلیسی؛
- `v4_validation/MARTINI3_PATH.md`: مراحل لازم برای مدل واقعی Martini 3.

## امانت‌داری علمی
روی GitHub باید مدل را **Exploratory reduced-unit casein-inspired coarse-grained MD** معرفی کرد. از عبارت‌های *Martini 3 MD performed*, *experimentally validated SAXS*, *calcium binding quantified* استفاده نکن.
