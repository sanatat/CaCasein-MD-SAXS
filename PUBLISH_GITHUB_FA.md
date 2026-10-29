# انتشار CaCasein-MD-SAXS روی GitHub

مخزن عمومی موجود: https://github.com/sanatat/CaCasein-MD-SAXS

این مخزن هنگام تهیه بسته **خالی** بود. پروژه را برای بارگذاری روی همین مخزن تنظیم کرده‌ایم؛ نیازی به ایجاد مخزن جدید نیست.

## روش ساده با Linux Mint

1. فایل ZIP نهایی را داخل `~/Downloads` استخراج کن.
2. اگر Git نصب نیست: `sudo apt update && sudo apt install git`.
3. هویت Git خودت (نه رمز عبور GitHub) را تنظیم کن، اگر قبلاً تنظیم نکرده‌ای:

```bash
git config --global user.name "Your Name"
git config --global user.email "your-email@example.com"
```

4. وارد پوشه پروژه شو و اسکریپت انتشار را اجرا کن:

```bash
cd ~/Downloads/CaCasein-MD-SAXS-GitHub-ready
bash PUBLISH_TO_GITHUB.sh
```

اسکریپت از Git استفاده می‌کند و برای مخزن موجود `origin` تنظیم می‌کند؛ **force push نمی‌کند**. GitHub دیگر پسورد عادی حساب را برای HTTPS push قبول نمی‌کند. برای ورود می‌توانی GitHub CLI نصب و اجرا کنی (`sudo apt install gh` سپس `gh auth login`) و بعد اسکریپت را اجرا کنی. توکن یا رمز عبور را در چت ارسال نکن.

بعد از push به https://github.com/sanatat/CaCasein-MD-SAXS برو و تب `Actions` را برای نتیجه تست‌ها بررسی کن.

## اجرای تست محلی قبل از انتشار

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python -m pytest -q tests v2_sequence_calcium/tests v3_explicit_ca/tests v4_validation/tests v5_saxs_identifiability/tests
```

**دقت علمی:** SAXS مصنوعی است؛ نسخه‌های V1–V4 مدل‌های ساده‌شده‌اند؛ V5 مدل جدید Martini 3 یا تأیید تجربی نیست. این پروژه به‌عنوان portfolio روش‌های محاسباتی منتشر می‌شود، نه یافته تأییدشده در زیست‌فیزیک کازئین.
