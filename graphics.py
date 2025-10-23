import pandas as pd
import matplotlib.pyplot as plt
import os
import re
import numpy as np
from typing import List, Dict

def plot_comparison_graphs(csv_dir: str = "./notes/res4/csvs", output_dir: str = "./plots"):
    """
    Строит графики сравнения тестов для фиксированных dist_name и типа fixed_value
    """
    
    # Создаем директорию для графиков
    os.makedirs(output_dir, exist_ok=True)
    
    # Получаем все CSV файлы
    csv_files = [f for f in os.listdir(csv_dir) if f.endswith('.csv')]
    
    # Группируем файлы по dist_name и типу fixed_value
    file_groups = {}
    
    for file in csv_files:
        # Парсим имя файла: {test_name}_{dist_name}_{h1|h2=value}
        match = re.match(r'([A-Z]+)_([A-Z]+)_(h[12])=([\d.]+)\.csv', file)
        if match:
            test_name, dist_name, fixed_type, fixed_value = match.groups()
            key = (dist_name, fixed_type)
            if key not in file_groups:
                file_groups[key] = []
            file_groups[key].append((test_name, float(fixed_value), file))
    
    # Строим графики для каждой группы
    for (dist_name, fixed_type), files in file_groups.items():
        print(f"Обрабатываю {dist_name} с фиксированным {fixed_type}...")
        
        # Читаем все DF для этой группы
        dfs = {}
        varying_param = 'h1' if fixed_type == 'h2' else 'h2'
        
        for test_name, fixed_value, filename in files:
            df = pd.read_csv(os.path.join(csv_dir, filename))
            # Обрабатываем названия колонок
            if fixed_type == 'h2':
                df.columns = [col.replace('h1...', f'{varying_param}=') for col in df.columns]
            else:
                df.columns = [col.replace('h2...', f'{varying_param}=') for col in df.columns]
            
            # Сохраняем первую колонку как индекс
            df = df.set_index(df.columns[0])
            dfs[(test_name, fixed_value)] = df

        
        new_out_dir = output_dir + f"/{dist_name}_fixed_{fixed_type}"
        os.makedirs(new_out_dir, exist_ok=1)
        os.makedirs(new_out_dir + '/fixed_n', exist_ok=1)
        os.makedirs(new_out_dir + f'/fixed_{varying_param}', exist_ok=1)
    
        # 3.1) Графики для строк (фиксируем n)
        plot_by_rows(dfs, dist_name, fixed_type, varying_param, new_out_dir + '/fixed_n')
        
        # 3.2) Графики для столбцов (фиксируем varying_param)
        plot_by_columns(dfs, dist_name, fixed_type, varying_param, new_out_dir + f'/fixed_{varying_param}')


def plot_by_rows(dfs: Dict, dist_name: str, fixed_type: str, varying_param: str, output_dir: str):
    """
    Строит графики для фиксированных значений n
    """
    
    # Получаем все уникальные значения n из индексов
    all_n_values = set()
    for df in dfs.values():
        all_n_values.update([idx for idx in df.index if '(n=' in idx])
    n_values = sorted(all_n_values, key=lambda x: int(x.split('(n=')[1].split(')')[0]))

    # Для каждого n строим отдельный график
    for n_value in n_values:
        plt.figure(figsize=(12, 8))
        
        # Собираем данные для каждого теста
        for (test_name, fixed_value), df in dfs.items():
            if n_value in df.index:
                row_data = df.loc[n_value]
                
                # Извлекаем значения varying_param и соответствующие значения мощности
                x_values = []
                y_values = []
                
                for col_name, value in row_data.items():
                    
                    if pd.notna(value) and f'{varying_param}=' in col_name:
                        param_value = float(col_name.split('=')[1])
                        x_values.append(param_value)
                        y_values.append(value)
                
                if x_values and y_values:
                    # Сортируем по x_values
                    sorted_data = sorted(zip(x_values, y_values))
                    x_sorted, y_sorted = zip(*sorted_data)
                    
                    plt.plot(x_sorted, y_sorted, 
                            marker='o', 
                            linewidth=2,
                            label=f'{test_name} ({fixed_type}={fixed_value})')
        
        plt.xlabel(f'{varying_param} values', fontsize=12)
        plt.ylabel('Empirical Power', fontsize=12)
        plt.title(f'{dist_name} Distribution with {fixed_type} = {fixed_value}\n  Fixed {n_value[n_value.index("(")+1:n_value.index(")")]}\n', fontsize=14)
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        
        # Сохраняем график
        filename = f"{dist_name}_{fixed_type}_n_{n_value.replace('=', '_')}.png"
        plt.savefig(os.path.join(output_dir, filename), dpi=300, bbox_inches='tight')
        plt.close()
        
        # print(f"Сохранен график: {filename}")

def plot_by_columns(dfs: Dict, dist_name: str, fixed_type: str, varying_param: str, output_dir: str):
    """
    Строит графики для фиксированных значений varying_param
    """
    
    # Получаем все уникальные значения varying_param из названий колонок
    all_param_values = set()
    for df in dfs.values():
        for col in df.columns:
            if f'{varying_param}=' in col:
                param_value = float(col.split('=')[1])
                all_param_values.add(param_value)
    
    param_values = sorted(all_param_values)


    # Для каждого значения varying_param строим отдельный график
    for param_value in param_values:
        plt.figure(figsize=(12, 8))
        
        column_name = f'{varying_param}={param_value}'
        
        # Собираем данные для каждого теста
        for (test_name, fixed_value), df in dfs.items():
            if column_name in df.columns:
                # Берем только строки с n значениями (исключаем asp_power)
                n_rows = [idx for idx in df.index if '(n=' in idx]
                n_data = df.loc[n_rows, column_name]

                if not n_data.empty:
                    # Извлекаем значения n и соответствующие значения мощности
                    n_values = [int(idx.split('(n=')[1].split(')')[0]) for idx in n_data.index]
                    powers = n_data.values
                    
                    # Сортируем по n_values
                    sorted_data = sorted(zip(n_values, powers))
                    n_sorted, powers_sorted = zip(*sorted_data)
                    
                    plt.plot(n_sorted, powers_sorted, 
                            marker='s', 
                            linewidth=2,
                            label=f'{test_name} ({fixed_type}={fixed_value})')
        
        plt.xlabel('Sample Size (n)', fontsize=12)
        plt.ylabel('Empirical Power', fontsize=12)
        plt.title(f'{dist_name} Distribution with {fixed_type} = {fixed_value} \n fixed {varying_param}={param_value}\n', fontsize=14)
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        
        # Сохраняем график
        filename = f"{dist_name}_{fixed_type}_{varying_param}_{param_value}.png"
        plt.savefig(os.path.join(output_dir, filename), dpi=300, bbox_inches='tight')
        plt.close()
        
        # print(f"Сохранен график: {filename}")


if __name__ == "__main__":

    plot_comparison_graphs(csv_dir = './notes/res4/csvs/', output_dir='./graphics/')
    

    
    print("Все графики построены и сохранены!")