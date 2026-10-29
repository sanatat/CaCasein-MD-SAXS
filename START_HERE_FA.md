# راهنمای اجرای پروژه روی لینوکس Mint / Ubuntu

این یک مدل فیزیکی کوچک از تجمع زنجیره‌های شبه‌پروتئینی است، نه شبیه‌سازی کازئین واقعی، یون کلسیم واقعی یا SAXS آزمایشگاهی.

```bash
unzip CaCasein-MD-SAXS.zip
cd CaCasein-MD-SAXS
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python -m pytest -q tests
python cacasein_md_saxs.py --steps 12000 --replicas 3 --stride 200 --out new_results
```

برای اجرای آزمایشی کوتاه:

```bash
python cacasein_md_saxs.py --steps 2000 --replicas 1 --stride 100 --out quick_result
```

برای دیدن نمودارها، پوشه `results/` یا `new_results/` را باز کنید؛ فایل `03_debye_saxs.png` منحنی پراکندگی مصنوعی و فایل `05_snapshot.png` ساختار نهایی زنجیره‌هاست.

نکته مهم: پارامتر `epsilon` فقط نیروی جاذبه مؤثر بین زنجیره‌هاست و معادل غلظت واقعی Ca2+ نیست. واحد طول و زمان هم بی‌بعد/کاهش‌یافته هستند؛ نباید آن‌ها را nm یا ns نام‌گذاری کرد.
