import re
import pandas as pd

def parse_file(filename):
    with open(filename, 'r', encoding='utf-8') as file:
        content = file.read()
    
    # Разделяем содержимое по h2
    sections = content.split('h2 = ')
    
    results = {}
    n_values = {}
    
    for section in sections[1:]:  # Пропускаем первую пустую секцию
        lines = section.strip().split('\n')
        
        # Извлекаем значение h2
        h2_value = float(lines[0])
        
        # Извлекаем значения n
        n_list = []
        emp_powers = None
        asp_power = None
        
        for line in lines:
            # Ищем строки с n
            if 'n = ' in line and 'Время выполнения:' in line:
                n_match = re.search(r'n = (\d+)', line)
                if n_match:
                    n_list.append(int(n_match.group(1)))
            elif line.startswith('emp_powers = '):
                # Извлекаем числа из emp_powers
                emp_match = re.findall(r'[\d.]+', line.split('=')[1])
                emp_powers = [float(x) for x in emp_match]
            elif line.startswith('asp_power = '):
                # Извлекаем asp_power
                asp_match = re.search(r'[\d.]+', line.split('=')[1])
                if asp_match:
                    asp_power = float(asp_match.group())
        
        if emp_powers is not None and n_list:
            # Сохраняем все значения
            all_values = emp_powers + [asp_power] if asp_power is not None else emp_powers
            results[h2_value] = all_values
            n_values[h2_value] = n_list
    
    return results, n_values

def create_transposed_dataframe_with_n(results, n_values):
    # Создаем списки для данных
    h2_values = sorted(results.keys())
    
    # Определяем количество emp_powers (все кроме последнего asp_power)
    num_emp_powers = len(results[h2_values[0]]) - 1
    
    # Проверяем, что количество n совпадает с количеством emp_powers
    for h2 in h2_values:
        if len(n_values[h2]) != num_emp_powers:
            print(f"Предупреждение: для h2={h2} количество n ({len(n_values[h2])}) не совпадает с количеством emp_powers ({num_emp_powers})")
    
    # Создаем имена строк с n
    index_names = []
    data = []
    
    # Добавляем emp_powers с n
    for i in range(num_emp_powers):
        row_name = f"EP (n={n_values[h2_values[0]][i]})"
        index_names.append(row_name)
        
        row_data = []
        for h2 in h2_values:
            if i < len(results[h2]) - 1:  # -1 потому что последний элемент это asp_power
                row_data.append(results[h2][i])
            else:
                row_data.append(None)
        data.append(row_data)
    
    # Добавляем asp_power как отдельную строку
    index_names.append('AP')
    asp_power_row = [results[h2][-1] for h2 in h2_values]  # последний элемент каждого списка
    data.append(asp_power_row)
    
    # Создаем DataFrame
    df = pd.DataFrame(data, index=index_names, columns=[f'h2 = {h2}' for h2 in h2_values])
    
    return df

def main():
    filename = './notes/cauchy_res.txt'  # Замените на имя вашего файла
    
  
    results, n_values = parse_file(filename)
    df = create_transposed_dataframe_with_n(results, n_values)
    df.to_csv('./tex_shit/df.csv', columns=df.columns)
    
    # f = open('./table.tex', 'w')
    # f.writelines([
    #     "\\documentclass{article}\n",
    #     "\\usepackage[utf8]{inputenc}\n",
    #     "\\usepackage[T2A]{fontenc}\n",
    #     "\\usepackage[russian]{babel}\n",
    #     "\\usepackage{booktabs}\n",
    #     "\\usepackage{amsmath}\n",
    #     "\\usepackage{graphicx}\n",
    #     "\\usepackage{float}\n",
    #     "\\begin{document}\n"
    # ])
    # df.to_latex(buf = f, float_format="%.3f", bold_rows=True, caption='ABOBA', position='h!')
    # f.write("\\end{document}")
    # f.close()

        


if __name__ == "__main__":
    main()