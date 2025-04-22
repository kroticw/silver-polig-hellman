#include <stdio.h>
#include <stdlib.h>
#include <math.h> // Для pow в factorize и sqrt

// Структура для хранения простого множителя и его степени
typedef struct {
    unsigned long long p; // Простой множитель
    unsigned int alpha;   // Степень
} PrimeFactor;

// Прототипы вспомогательных функций
unsigned long long power(unsigned long long base, unsigned long long exp, unsigned long long modulus);
int factorize(unsigned long long n_orig, PrimeFactor** factors_ptr, int* count_ptr);
unsigned long long discrete_log_prime_power(
    unsigned long long a,
    unsigned long long b,
    unsigned long long p,
    unsigned int alpha,
    unsigned long long q
    );
unsigned long long solve_crt(const unsigned long long* remainders, const unsigned long long* moduli, int n);
unsigned long long modInverse(unsigned long long n, unsigned long long modulus); // Для CRT

// Функция модульного возведения в степень (base^exp) % modulus
// Использует алгоритм бинарного возведения в степень
unsigned long long power(unsigned long long base, unsigned long long exp, unsigned long long modulus) {
    unsigned long long res = 1;
    base %= modulus;
    while (exp > 0) {
        if (exp % 2 == 1) res = (__int128)res * base % modulus; // Используем __int128 для промежуточного результата во избежание переполнения
        base = (__int128)base * base % modulus;
        exp /= 2;
    }
    return res;
}

// Функция для разложения числа n на простые множители
// Возвращает 0 при успехе, -1 при ошибке (например, не удалось выделить память)
// factors_ptr - указатель на массив PrimeFactor, count_ptr - количество найденных факторов
// Простая реализация пробным делением.
int factorize(const unsigned long long n_orig, PrimeFactor** factors_ptr, int* count_ptr) {
    *count_ptr = 0;
    *factors_ptr = NULL;
    int capacity = 10; // Начальная вместимость
    PrimeFactor* temp_factors = malloc(capacity * sizeof(PrimeFactor));
    if (!temp_factors) return -1; // Ошибка выделения памяти

    unsigned long long temp_n = n_orig;

    // Обработка множителя 2
    if (temp_n % 2 == 0) {
        if (*count_ptr >= capacity) {
             capacity *= 2;
             PrimeFactor* new_factors = (PrimeFactor*)realloc(temp_factors, capacity * sizeof(PrimeFactor));
             if (!new_factors) { free(temp_factors); return -1; }
             temp_factors = new_factors;
        }
        temp_factors[*count_ptr].p = 2;
        temp_factors[*count_ptr].alpha = 0;
        while (temp_n % 2 == 0) {
            temp_factors[*count_ptr].alpha++;
            temp_n /= 2;
        }
        (*count_ptr)++;
    }

    // Обработка нечетных множителей
    // Делим до sqrt(temp_n)
    unsigned long long limit = 1;
    // Используем __int128 для безопасного вычисления limit*limit
    while ((__int128)limit * limit <= temp_n) {
        limit++;
    }
    limit--; // вернулись к значению, чей квадрат <= temp_n

    for (unsigned long long i = 3; i <= limit; i += 2) {
        if (temp_n % i == 0) {
             if (*count_ptr >= capacity) {
                 capacity *= 2;
                 PrimeFactor* new_factors = (PrimeFactor*)realloc(temp_factors, capacity * sizeof(PrimeFactor));
                 if (!new_factors) { free(temp_factors); return -1; }
                 temp_factors = new_factors;
             }
            temp_factors[*count_ptr].p = i;
            temp_factors[*count_ptr].alpha = 0;
            while (temp_n % i == 0) {
                temp_factors[*count_ptr].alpha++;
                temp_n /= i;
            }
            (*count_ptr)++;
            // Обновляем предел после деления
            limit = 1;
             while ((__int128)limit * limit <= temp_n) {
                 limit++;
             }
             limit--;
        }
    }

    // Если temp_n стал простым числом > sqrt(исходного n)
    if (temp_n > 1) {
        if (*count_ptr >= capacity) {
             capacity++; // Достаточно увеличить на 1
             PrimeFactor* new_factors = (PrimeFactor*)realloc(temp_factors, capacity * sizeof(PrimeFactor));
             if (!new_factors) { free(temp_factors); return -1; }
             temp_factors = new_factors;
        }
        temp_factors[*count_ptr].p = temp_n;
        temp_factors[*count_ptr].alpha = 1;
        (*count_ptr)++;
    }

    *factors_ptr = temp_factors;
    return 0; // Успех
}

// Расширенный алгоритм Евклида для нахождения НОД(a, b) и коэффициентов x, y
// таких что ax + by = НОД(a, b)
// Возвращает НОД(a, b)
// *x и *y будут содержать коэффициенты Безу
// используем __int128 для избежания переполнения при промежуточных вычислениях
long long extended_gcd(const long long a, const long long b, long long *x, long long *y) {
    if (a == 0) {
        *x = 0;
        *y = 1;
        return b;
    }

    long long x1, y1;
    const long long gcd = extended_gcd(b % a, a, &x1, &y1);

    *x = y1 - b / a * x1; // Приведение типа для деления
    *y = x1;

    return gcd;
}

// Функция для нахождения модульного обратного n по модулю modulus
// использует расширенный алгоритм Евклида,
// возвращает n^{-1} mod modulus, или 0, если обратного не существует (n и modulus не взаимно просты)
unsigned long long modInverse(const unsigned long long n, const unsigned long long modulus) {
    long long x, y;
    const long long n_ll = (long long)n;
    const long long mod_ll = (long long)modulus;

    const long long g = extended_gcd(n_ll, mod_ll, &x, &y);
    if (g != 1) {
        // Обратный элемент не существует
        fprintf(stderr, "Ошибка: Обратный элемент для %llu по модулю %llu не существует.", n, modulus);
        return 0; // Или можно использовать другое значение для индикации ошибки
    }
    // Делаем x положительным
    const unsigned long long res = (x % mod_ll + mod_ll) % mod_ll;
    return res;
}

// Функция для вычисления дискретного логарифма x = log_a(b) mod p^alpha
// где p - простой множитель q-1
unsigned long long discrete_log_prime_power(
    const unsigned long long a,
    const unsigned long long b,
    const unsigned long long p,
    const unsigned int alpha,
    const unsigned long long q
    )
{
    unsigned long long p_pow_alpha = 1;
    for(int i = 0; i < alpha; ++i) p_pow_alpha *= p;

    unsigned long long x = 0;
    unsigned long long p_pow_k = 1; // p^k

    // Предварительно вычислим a_inv = a^{-1} mod q
    const unsigned long long a_inv = modInverse(a, q);
    if (a_inv == 0) { // Проверка на случай, если обратный элемент не найден
        fprintf(stderr, "Ошибка: не удалось найти обратный элемент для a=%llu mod q=%llu\n", a, q);
        return -1;
    }

    // Вычисляем x = x_0 + x_1*p + ... + x_{alpha-1}*p^{alpha-1}
    unsigned long long current_b = b; // b_k в обозначениях из описания
    for (unsigned int k = 0; k < alpha; ++k) {
        // Вычисляем h = a^((q-1)/p)
        const unsigned long long h = power(a, (q - 1) / p, q);

        // Вычисляем c = b_k^((q-1)/p^{k+1}) mod q
        const unsigned long long p_pow_k_plus_1 = p_pow_k * p;
        const unsigned long long exponent = (q - 1) / p_pow_k_plus_1; // степень -- ((q-1)/p^{k+1})
        const unsigned long long c = power(current_b, exponent, q);

        // Находим x_k путем перебора j от 0 до p-1
        // Ищем j такое, что h^j = c mod q
        unsigned long long x_k = -1; // Значение по умолчанию, если не найдено
        unsigned long long h_pow_j = 1; // h^j
        for (unsigned long long j = 0; j < p; ++j) {
            if (h_pow_j == c) {
                x_k = j;
                break;
            }
            h_pow_j = (__int128)h_pow_j * h % q;
        }

        if (x_k == (unsigned long long)-1) {
             fprintf(stderr, "Ошибка: не удалось найти x_%u для p=%llu\n", k, p);
             // Это может произойти, если b не является степенью a по модулю q,
             // или из-за ограничений unsigned long long
             return -1; // Индикация ошибки
        }

        // Добавляем найденный коэффициент к результату
        x = (x + (__int128)x_k * p_pow_k) % p_pow_alpha; // Собираем результат по модулю p^alpha

        // Обновляем b_k для следующей итерации: b_{k+1} = b_k * (a^{-1})^{x_k * p^k}
        // b_{k+1} = b_k * (a^{-x_k})^{p^k}
        const unsigned long long correction_exp = (__int128)x_k * p_pow_k % (q-1); // Вычисляем показатель степени для коррекции по модулю порядка группы (q-1)
        const unsigned long long correction_term = power(a_inv, correction_exp, q);
        current_b = (__int128)current_b * correction_term % q;

        // Обновляем p^k
        p_pow_k *= p;

    }

    return x;
}

// Функция для решения системы сравнений с помощью Китайской теоремы об остатках (CRT)
// remainders[i] = x mod moduli[i]
// moduli[i] должны быть попарно взаимно простыми (что верно для p_i^alpha_i)
unsigned long long solve_crt(const unsigned long long* remainders, const unsigned long long* moduli, const int n) {
    unsigned long long M = 1;
    for (int i = 0; i < n; i++) {
        // Проверка на переполнение при вычислении M
        // используем __int128 для промежуточного хранения произведения
        if ((__int128)M * moduli[i] > (unsigned long long)-1) {
             fprintf(stderr, "Предупреждение: Возможно переполнение при вычислении общего модуля M в CRT.\n");
             // Продолжаем, но результат может быть неверным
        }
        M *= moduli[i];
    }

    unsigned long long result = 0;

    for (int i = 0; i < n; i++) {
        const unsigned long long a_i = remainders[i];
        const unsigned long long m_i = moduli[i];
        const unsigned long long M_i = M / m_i;
        const unsigned long long y_i = modInverse(M_i, m_i);

        if (y_i == 0) {
             fprintf(stderr, "Ошибка: не удалось найти обратный элемент в CRT для M_i=%llu mod m_i=%llu\n", M_i, m_i);
             return -1; // Ошибка CRT
        }

        // Вычисляем i-й член суммы с использованием __int128 для предотвращения переполнения
        __int128 term = (__int128)a_i * M_i % M; // Сначала по модулю M, чтобы уменьшить числа
        term = term * y_i % M;

        result = (result + term) % M;
    }

    return result;
}

int main(const int argc, char *argv[]) {
    // x = log_a_b mod q

    // x = log_2_62 mod 181
    const unsigned long long a = 2;
    const unsigned long long b = 62;
    const unsigned long long q = 181;

    // x = log_6_7531 mod 8101
    // const unsigned long long a = 6;
    // const unsigned long long b = 7531;
    // const unsigned long long q = 8101;

    // x = log_2_28 mod 37 -- пример с википедии
    // const unsigned long long a = 2;
    // const unsigned long long b = 28;
    // const unsigned long long q = 37;

    printf("Вычисление x = log_%llu(%llu) mod %llu\n", a, b, q);

    // Шаг 1: Разложить q-1 на простые множители
    const unsigned long long q_minus_1 = q - 1;
    PrimeFactor* factors = NULL;
    int factor_count = 0;
    printf("Факторизация q-1 = %llu...\n", q_minus_1);
    if (factorize(q_minus_1, &factors, &factor_count) != 0 || factor_count == 0) {
        fprintf(stderr, "Ошибка: Не удалось факторизовать q-1.\n");
        free(factors); // Освобождаем память, даже если она NULL
        return 1;
    }

    printf("Найдено %d простых множителей q-1:\n", factor_count);
    for (int i = 0; i < factor_count; i++) {
        printf("  p_%d = %llu, alpha_%d = %u\n", i, factors[i].p, i, factors[i].alpha);
    }

    // Шаг 2 и 3: Вычислить x mod p_i^alpha_i для каждого множителя и собрать результаты
    unsigned long long* r = malloc(factor_count * sizeof(unsigned long long));
    unsigned long long* moduli = malloc(factor_count * sizeof(unsigned long long));
    if (!r || !moduli) {
        fprintf(stderr, "Ошибка: Не удалось выделить память для результатов CRT.\n");
        free(factors);
        free(r);
        free(moduli);
        return 1;
    }

    printf("Вычисление дискретных логарифмов по модулю степеней простых чисел...\n");
    for (int i = 0; i < factor_count; i++) {
        const unsigned long long p_i = factors[i].p;
        const unsigned int alpha_i = factors[i].alpha;
        unsigned long long p_pow_alpha = 1;
        for(int j=0; j<alpha_i; ++j) p_pow_alpha *= p_i;

        printf("  Вычисление x mod %llu^%u = %llu...\n", p_i, alpha_i, p_pow_alpha);
        const unsigned long long x_i = discrete_log_prime_power(a, b, p_i, alpha_i, q);

        if (x_i == (unsigned long long)-1) {
            fprintf(stderr, "Ошибка при вычислении дискретного логарифма для p=%llu.\n", p_i);
            free(factors);
            free(r);
            free(moduli);
            return 1;
        }

        r[i] = x_i;
        moduli[i] = p_pow_alpha;
        printf("    x === %llu (mod %llu)\n", r[i], moduli[i]);
    }

    // Шаг 3-2: Использовать Китайскую теорему об остатках (CRT)
    printf("Применение Китайской теоремы об остатках...\n");
    const unsigned long long result_x = solve_crt(r, moduli, factor_count);

    if (result_x == (unsigned long long)-1) {
         fprintf(stderr, "Ошибка при применении CRT.\n");
         free(factors);
         free(r);
         free(moduli);
         return 1;
    }

    printf("Проверка: %llu^%llu mod %llu =? %llu\n", a, result_x, q, b);
    const unsigned long long check = power(a, result_x, q);
    if (check == b) {
        printf("Результат: x = %llu\n", result_x);
    } else {
        printf("Проверка не удалась. Вычислено: %llu. Возможно, дискретный логарифм не существует или произошла ошибка вычислений.\n", check);
    }

    // Очистка памяти
    free(factors);
    free(r);
    free(moduli);

    return 0;
}
