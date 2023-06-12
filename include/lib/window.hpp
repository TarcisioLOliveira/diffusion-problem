/*
 *   Copyright (C) 2023 Tarcísio Ladeia de Oliveira.
 *
 *   This file is part of diffusion-problem
 *
 *   diffusion-problem is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   diffusion-problem is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with diffusion-problem.  If not, see <https://www.gnu.org/licenses/>.
 *
 */

#ifndef DPLIB_WINDOW_HPP
#define DPLIB_WINDOW_HPP

#include <SFML/Config.hpp>
#include <SFML/Graphics/Color.hpp>
#include <SFML/Graphics/Sprite.hpp>
#include <SFML/Graphics/Texture.hpp>
#include <SFML/System/Vector2.hpp>
#include <SFML/Window.hpp>
#include <SFML/Graphics.hpp>

namespace dplib{

class Window{
    public:
    enum class ColorMap{
        GRAYSCALE,
        HSV
    };
    Window(size_t window_width, size_t window_height, size_t mesh_width, size_t mesh_height, std::string name);

    void update();
    void update(const std::vector<double>& mesh);
    void update(const std::vector<double>& mesh, const double min_x, const double max_x);
    void save_image();

    inline bool is_open(){
        return this->window.isOpen();
    }

    private:
    size_t window_width, window_height, W, H, legend_width, legend_height;
    sf::Texture img;
    sf::Texture legend_img;
    sf::Sprite sprite;
    sf::Sprite legend;
    sf::RenderWindow window;
    std::vector<sf::Uint8> pixels;
    std::string name;
    sf::Text text_max;
    sf::Text text_min;
    sf::Font font;
    ColorMap colormap = ColorMap::GRAYSCALE;

    void to_grayscale(const std::vector<double>& mesh, const double min_x, const double max_x);
    void to_hsv(const std::vector<double>& mesh, const double min_x, const double max_x);
    void update(const double min_x, const double max_x);

    void legend_grayscale();
};

}

#endif
